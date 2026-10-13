"""绘图专用模块锁、子进程退出和可核验临时目录收尾。"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time
from datetime import datetime, timezone


def _sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _files(root, scientific=False):
    if not root.is_dir():
        return {}
    return {str(p.relative_to(root)): _sha(p) for p in sorted(root.rglob('*'))
            if p.is_file() and (not scientific or
            ('00_audit' not in p.relative_to(root).parts and '01_audit' not in p.relative_to(root).parts
             and 'logs' not in p.relative_to(root).parts
             and p.suffix.lower() not in {'.pdf', '.png', '.html', '.md', '.log'}
             and not p.name.endswith('_plot_source.tsv')))}


class PlotRun:
    """同一分析模块只有一个绘图进程；子程序继承 flock 文件描述符。"""

    def __init__(self, root, evidence_dir=None, result_name='result'):
        self.root = Path(root).resolve()
        if Path(result_name).name != result_name:
            raise ValueError('结果目录必须是模块根的直接子目录名')
        self.result_name = result_name
        if evidence_dir is None:
            evidence_dir = self.root / 'logs/plot_failures'
        self.evidence_dir = Path(evidence_dir)
        self.children = []
        self.work = []
        self.handlers = {}
        self.logs = []
        self.commands = []

    def __enter__(self):
        self.lock = (self.root / '.plot.lock').open('a+')
        try:
            fcntl.flock(self.lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BaseException:
            self.lock.close()
            raise RuntimeError(f'已有绘图父进程或子进程持有模块锁: {self.root}') from None
        try:
            self._recover()
            for sig in (signal.SIGTERM, signal.SIGINT):
                self.handlers[sig] = signal.signal(sig, self._interrupt)
            return self
        except BaseException:
            self.lock.close()
            raise

    @staticmethod
    def _interrupt(signum, frame):
        raise InterruptedError(f'绘图收到终止信号 {signum}')

    def stage(self, prefix):
        marker = {'version': 1, 'pid': os.getpid(), 'result': self.result_name,
                  'source_science': _files(self.root / self.result_name, scientific=True)}
        oldmask = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGTERM, signal.SIGINT})
        try:
            work = Path(tempfile.mkdtemp(prefix=prefix, dir=self.root))
            self.work.append(work)
            (work / '.plot-stage.json').write_text(json.dumps(marker, ensure_ascii=False, indent=2) + '\n')
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, oldmask)
        return work

    def command(self, argv, *, check=False, **kwargs):
        """等待独立进程组；TERM/INT 先停止并 wait，锁不会提前关闭。"""
        kwargs.pop('start_new_session', None)
        inherited = tuple(kwargs.pop('pass_fds', ())) + (self.lock.fileno(),)
        if kwargs.get('stdout') is not None and hasattr(kwargs['stdout'], 'name'):
            path = Path(kwargs['stdout'].name)
            if path.is_file():
                self.logs.append(path)
        if 'stdout' not in kwargs and self.work:
            log = self.work[-1] / f'command_{len(self.commands) + 1:02d}.log'
            with log.open('wb') as handle:
                self.logs.append(log)
                kwargs.setdefault('stderr', subprocess.STDOUT)
                return self.command(argv, check=check, stdout=handle, **kwargs)
        child = subprocess.Popen(argv, start_new_session=True, pass_fds=inherited, **kwargs)
        self.children.append(child)
        record = {'argv': [str(value) for value in argv], 'pid': child.pid, 'returncode': None}
        self.commands.append(record)
        try:
            result = subprocess.CompletedProcess(argv, child.wait())
            record['returncode'] = result.returncode
        except BaseException:
            self._stop_children()
            raise
        if check:
            result.check_returncode()
        return result

    def _stop_children(self):
        # Children can have grandchildren in the same session. Signal each group
        # even when its direct leader has already exited.
        for child in self.children:
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        for child in self.children:
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        # A leader may exit before a TERM-ignoring grandchild; ensure the group
        # cannot continue to write after cleanup or after releasing the lock.
        for child in self.children:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        for _ in range(100):
            live = []
            for process in Path('/proc').iterdir():
                if not process.name.isdigit():
                    continue
                try:
                    fields = (process / 'stat').read_text().rsplit(')', 1)[1].split()
                    if fields[0] != 'Z' and int(fields[2]) in {child.pid for child in self.children}:
                        live.append(process.name)
                except (OSError, ValueError):
                    pass
            if not live:
                return
            time.sleep(.05)
        raise RuntimeError(f'绘图子进程组仍未退出，禁止清理: {live}')

    def publish(self, stage, result):
        stage, result = Path(stage), Path(result)
        if result != self.root / self.result_name:
            raise ValueError('发布目标不是本次模块锁所属结果目录')
        work = next(p for p in self.work if stage.is_relative_to(p))
        marker_path = work / '.plot-stage.json'
        marker = json.loads(marker_path.read_text())
        marker['publication_expected'] = _files(stage)
        marker_path.write_text(json.dumps(marker, ensure_ascii=False, indent=2) + '\n')
        previous = work / 'previous_result'
        # Block catchable signals across the two renames. SIGKILL recovery uses
        # previous_result plus the exact complete output identity above.
        oldmask = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGTERM, signal.SIGINT})
        try:
            if result.exists():
                result.rename(previous)
            try:
                stage.rename(result)
            except BaseException:
                if previous.exists() and not result.exists():
                    previous.rename(result)
                raise
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, oldmask)
        if previous.exists():
            shutil.rmtree(previous)

    def _remove(self, work):
        marker_path = work / '.plot-stage.json'
        marker = json.loads(marker_path.read_text())
        result = self.root / marker['result']
        previous = work / 'previous_result'
        if previous.exists():
            if not result.exists():
                previous.rename(result)
            elif _files(result) == marker.get('publication_expected'):
                shutil.rmtree(previous)
            else:
                raise RuntimeError(f'保留唯一旧结果：当前发布身份不完整: {work}')
        if marker.get('publication_expected') == _files(result):
            shutil.rmtree(work)
            return
        if _files(result, scientific=True) != marker['source_science']:
            raise RuntimeError(f'正式科学源身份改变，保留临时科学副本待核对: {work}')
        shutil.rmtree(work)

    def _recover(self):
        for work in sorted(self.root.iterdir()):
            if work.is_dir() and (work / '.plot-stage.json').is_file():
                for process in Path('/proc').iterdir():
                    if not process.name.isdigit() or int(process.name) == os.getpid():
                        continue
                    try:
                        command = (process / 'cmdline').read_bytes()
                        cwd = (process / 'cwd').resolve()
                    except (OSError, RuntimeError):
                        continue
                    if str(work).encode() in command or cwd.is_relative_to(work):
                        raise RuntimeError(f'临时目录仍被进程引用，禁止回收: {work}, PID {process.name}')
                self._remove(work)

    def _failure(self, error):
        target = self.evidence_dir / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ') + f'-{os.getpid()}')
        target.mkdir(parents=True, exist_ok=False)
        for index, log in enumerate(dict.fromkeys(self.logs)):
            if log.is_file():
                # Keep the relevant tail, never an entire scientific result copy.
                with log.open('rb') as handle:
                    handle.seek(max(0, log.stat().st_size - 8 * 1024 * 1024))
                    (target / f'{index:02d}_{log.name}').write_bytes(handle.read())
        for index, work in enumerate(self.work):
            marker = work / '.plot-stage.json'
            if marker.is_file():
                shutil.copy2(marker, target / f'{index:02d}_stage_signature.json')
        (target / 'failure.json').write_text(json.dumps({'status': 'FAILED', 'reason': str(error),
            'module': str(self.root), 'commands': self.commands,
            'stages': [str(p) for p in self.work]}, ensure_ascii=False, indent=2) + '\n')

    def __exit__(self, kind, error, traceback):
        # Ignore additional signals during deterministic child and stage cleanup.
        for sig in self.handlers:
            signal.signal(sig, signal.SIG_IGN)
        cleanup_error = None
        try:
            self._stop_children()
            if error is not None:
                self._failure(error)
            for work in self.work:
                if work.exists():
                    try:
                        self._remove(work)
                    except BaseException as issue:
                        cleanup_error = issue
            if cleanup_error is not None and error is None:
                raise cleanup_error
        finally:
            for sig, handler in self.handlers.items():
                signal.signal(sig, handler)
            self.lock.close()
        return False
