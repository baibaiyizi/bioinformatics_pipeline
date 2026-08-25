# scATAC / single细胞 multiome 最小闭环

该入口从已经生成的 RNA 与 ATAC count-level `.h5ad` 开始，完成模态分离 QC、共享细胞核对、ATAC TF-IDF/LSI 和可选 MuData 导出。原始 FASTQ/fragments 到 peak matrix 的过程仍应由 Cell Ranger ARC、ArchR、SnapATAC2 等成熟流程完成；本目录不复制这些仓库。

```bash
python scripts/run_multiome.py --config config/example.yml
python -m unittest discover -s tests -v
```

硬边界：

- RNA 与 ATAC 必须保留整数型原始 counts，且 cell barcode 唯一。
- paired 模式要求两模态细胞完全匹配；不静默取交集。
- `sample_column` 必须存在；下游 DE/DA 以样本 pseudobulk 为默认推断单位。
- TF-IDF、LSI 和联合表示是细胞级表征，不产生额外生物学重复。
- MuData/MultiVI 是可选后续层；未安装 `mudata/scvi-tools` 或没有 GPU 时会明确记录，不伪装成已经完成整合。
