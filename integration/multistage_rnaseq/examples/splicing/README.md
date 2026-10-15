# 剪接轨迹合成示例

运行 `bash run_example.sh`。S1故意把rMATS的b1设为MY、S2设为MO，用于验证统一的 `MO - MY` 方向；两个A3SS事件共享旧式部分坐标，但flanking坐标不同，用于验证完整坐标键不会碰撞。示例还保留一个PSI `NA`，并故意不提供BAM，预期Sashimi状态为 `missing_bam` 而非伪造图。

