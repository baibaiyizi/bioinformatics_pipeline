"""共同成图样式。项目运行自己的冻结副本，不依赖中心模板运行时。"""
from pathlib import Path
import json

WIDTH_MM = 183
DPI = 300
FONT = "Liberation Sans"
TEXT_PT = 8
TICK_PT = 7
PROFILE_RELATIVE = "../config/figure_style.json"
COLORS = {}
DIVERGING = []
PROFILE = {}


def load_profile(profile_path):
    path = Path(profile_path).resolve(strict=True)
    profile = json.loads(path.read_text())
    required = {"group", "de", "neutral", "diverging", "sequential", "binary", "quadrant"}
    if not required <= profile["palettes"].keys():
        raise ValueError("绘图profile缺少必需语义调色板")
    return profile


def configure_profile(profile_path):
    global WIDTH_MM, DPI, FONT, TEXT_PT, TICK_PT
    profile = load_profile(profile_path)
    PROFILE.clear()
    PROFILE.update(profile)
    COLORS.clear()
    COLORS.update(profile["palettes"]["group"])
    COLORS["NS"] = profile["palettes"]["de"]["ns"]
    DIVERGING[:] = list(profile["palettes"]["diverging"].values())
    typography = profile["typography"]
    WIDTH_MM, DPI, FONT = typography["width_mm"], typography["dpi"], typography["family"]
    TEXT_PT, TICK_PT = typography["text_pt"], typography["tick_pt"]
    return profile


configure_profile(Path(__file__).resolve().parent / PROFILE_RELATIVE)


def apply_style():
    import matplotlib as mpl
    from matplotlib import font_manager
    font_manager.findfont(FONT, fallback_to_default=False)
    mpl.rcParams.update({
        "font.family": FONT, "font.size": TEXT_PT, "axes.titlesize": TEXT_PT,
        "axes.labelsize": TEXT_PT, "xtick.labelsize": TICK_PT,
        "ytick.labelsize": TICK_PT, "legend.fontsize": TICK_PT,
        "legend.title_fontsize": TEXT_PT, "axes.grid": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": PROFILE["typography"]["line_mm"] * 72 / 25.4, "pdf.fonttype": 42,
        "ps.fonttype": 42, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.transparent": False,
        "savefig.dpi": DPI,
    })


def save_figure(fig, pdf, height_mm=None):
    import matplotlib.text
    from PIL import Image
    apply_style()
    pdf = Path(pdf)
    if pdf.suffix != ".pdf":
        raise ValueError("矢量主文件必须使用 .pdf")
    pdf.parent.mkdir(parents=True, exist_ok=True)
    if height_mm is None:
        height_mm = WIDTH_MM * fig.get_figheight() / fig.get_figwidth()
    fig.set_size_inches(WIDTH_MM / 25.4, height_mm / 25.4, forward=True)
    for item in fig.findobj(match=matplotlib.text.Text):
        item.set_fontfamily(FONT)
        item.set_fontsize(max(TICK_PT, min(TEXT_PT, item.get_fontsize())))
    for ax in fig.axes:
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.savefig(pdf, facecolor="white", transparent=False)
    png = pdf.with_suffix(".png")
    fig.savefig(png, dpi=DPI, facecolor="white", transparent=False)
    with Image.open(png) as source:
        rgb = source.convert("RGB")
        rgb.save(png, dpi=(DPI, DPI))
    return {"pdf": str(pdf), "png": str(png)}
