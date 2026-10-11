# 本地冻结实现

项目入口为上一层 `run.sh`。共享实现仅在中心 `templates/smallrnaseq` 维护；复制到项目后以源码 SHA256 冻结，不自动跟随模板更新。

物种、模式、设计与资源由所选数据集 `config/project.sh` 显式传入。`R/run_analysis.R` 只消费已解析的根与资源，不根据 cwd 或旧 PANDORA 布局猜测输入。SPORTS 上游保留原有物种剪切策略；R/Rmd 分析算法和分类边界不因目录调整而改变。

`preflight.py` 只读检查输入与模块适用性。完整共同分类、reads 守恒和统计资格检查仍在实际 R 执行中完成。`validate_figures.py` 检查生成文件的技术规格；不替代科学语义审阅。
