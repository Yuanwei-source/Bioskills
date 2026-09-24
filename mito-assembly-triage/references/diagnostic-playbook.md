# V2 动态诊断手册

诊断遵循 `INTAKE → HYPOTHESIZE → CHOOSE_TEST → EXECUTE → UPDATE → DECIDE → VERIFY → LEARN`。

| 观察 | 至少保留的替代解释 | 首个低成本检查 | 需要升级的证据 |
|---|---|---|---|
| 基因缺失/伪基因样命中 | 注释边界、真实缺失、contig 断裂/重复 | `seq_stats.py`、`annot_check.py`、独立同源命中 | reads/GFA 或独立注释 |
| CDS 移码/内部 stop | 边界/遗传密码表、碱基错误、真实伪基因 | 核对 table、坐标、翻译和质量 | pileup、独立蛋白/注释证据 |
| 接缝/控制区异常 | 重复错接、真实长度异质性、未闭合 | 端部重叠和局部比对 | 所有接缝的 reads/pair/GFA |
| 异常覆盖 | NUMT、污染、多倍型、组装错误 | contig/参考竞争比对和 MAPQ | 核参考、原始 reads、图 |
| 基因重排/顺序不同 | 真实重排、方向错误、参考不近缘 | 多参考/坐标和链方向检查 | 跨边界 reads 或长读长 |

只有检查能改变假设排序时才执行下一步。无法区分时停止并写出 `UNRESOLVED`。
