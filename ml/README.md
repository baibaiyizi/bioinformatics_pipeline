# SCI 通用预测建模最小闭环

这是一个面向**已经冻结且可追溯的样本级特征表**的通用分类/回归入口。它不负责从原始组学数据筛特征；特征冻结应由对应模态流程完成，并在 `feature_columns` 中显式记录。

核心保护：

- 缺失值填补、标准化和数据驱动特征选择全部位于 scikit-learn `Pipeline` 内。
- 有个体、窝、中心、批次等分组时，必须填写 `group_column`，外层和内层重采样均不拆分同一组。
- 调参与性能估计使用嵌套交叉验证；最终模型只在完成性能估计后于全部开发数据重拟合。
- 外部数据只做一次最终评估，不参与预处理拟合、特征选择或调参。
- 无独立外部数据时，报告会明确写为“内部验证”，不使用“已验证模型”。

## 运行

```bash
cd /home/h1028/workspace/model/ml
python scripts/run_ml.py --config config/example.yml
```

配置中的相对路径按配置文件所在目录解析。最小字段：

```yaml
input_csv: ../data/features.csv
target_column: outcome
task: classification
feature_columns: [gene_a, gene_b, age]
group_column: subject_id
output_dir: ../result/example
```

输出包括输入/环境审计、外层折叠预测和指标、校准表、决策曲线、置换重要性、最终模型、模型卡与可重跑配置。SHAP 是可选增强；未安装 `shap` 时会在模型卡中标记为未执行，不会用其他重要性冒充 SHAP。

## 输入契约

- 每行是一个统计单位；`id_column` 若提供必须唯一。
- `feature_columns` 必须显式给出且全部为数值列。
- 分类任务当前要求二分类；目标不得缺失。
- 分组交叉验证要求每个 fold 都具备可计算目标的组结构；样本过少时直接失败。
- 生存结局和深度学习不由该最小入口假装支持，应由 `sci_ml` 选择专用实现。

## 验证

```bash
python -m unittest discover -s tests -v
```

临床预测研究的写作按 TRIPOD+AI，方法与偏倚审查按 PROBAST+AI；这两个规范不替代数据集设计、校准和外部评估。
