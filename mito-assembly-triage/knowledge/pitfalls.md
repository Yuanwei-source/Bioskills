# 常见工具坑

| 工具 | 常见问题 | 处理 |
|------|----------|------|
| MITOS2 | 参考目录或依赖命令不在 PATH | 显式设置 `MITOS2_REFDIR`、`MITOS2_EXTRA_PATH` 并运行环境检查 |
| GetOrganelle | 重复区导致图复杂或无法环化 | 保留 scaffold，结合 reads 接缝证据选择修复路线 |
| MitoFinder | 参考格式、配置文件或 PATH 不完整 | 使用 GenBank 参考并检查工具环境 |
| BLAST | 短命中不代表完整基因边界 | 用全长参考、结构工具和翻译验证 |
| 坐标后处理 | 模糊碱基或边界被静默改写 | 每次替换都保留 reads 覆盖和工具证据 |
