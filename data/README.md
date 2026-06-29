# 数据准备说明

## 目录结构要求

```
Pan-Genome/                        # 你的泛基因组根目录
├── Fragaria_vesca_v1.0/
│   ├── Fvesca.ICL.txt             # HMMER hmmsearch --tblout 输出
│   ├── Fvesca.pep.fasta           # 蛋白组（或 *.fa / *.protein.fasta）
│   ├── genome.fa                  # 基因组序列（CDS 提取用）
│   └── annotation.gff             # 注释文件（CDS 提取用）
├── Fragaria_iinumae_v1.0/
│   └── ...
└── ...
```

## 关键注意事项

1. **HMMER 输出格式**: 必须用 `hmmsearch --tblout` 生成，默认读第 1 列 (ID) 和第 5 列 (E-value)
2. **FASTA 头匹配**: HMMER 结果中的 ID 必须能在蛋白 FASTA 的 `>` 头第一个空格分隔词匹配
3. **拟南芥参考**: 下载 TAIR10 蛋白组，ICL 参考基因为 `AT3G21720`
4. **路径名称**: 避免空格和中文特殊字符，推荐用下划线替代
