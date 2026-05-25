# Bioskills

这是一个面向生物信息学分析的智能 Skill 库，目标是把可复用的分析流程、脚本和 AI SOP (`SKILL.md`) 组织成清晰的技能模块，让 AI 助手能够更稳定地完成生信任务。

## 已集成 Skills

### 1. Enrichment
- **目录**: [Enrichment](./Enrichment/)
- **说明书**: [SKILL.md](./Enrichment/SKILL.md)
- **核心功能**:
  - 支持 ORA 和 GSEA 双模式。
  - 自动识别并转换基因 ID（Symbol、Ensembl、Entrez）。
  - 多物种支持（Human、Mouse、Pig、Rat、Cattle 等）。
  - 生成散点图、山脊图、Barcode 图等常见富集可视化结果。

### 2. BioFinder
- **目录**: [BioFinder](./BioFinder/)
- **说明书**: [SKILL.md](./BioFinder/SKILL.md)
- **核心功能**:
  - GEO-first 的公开组学数据发现与评估 workflow。
  - 支持候选检索、accession 级 evidence extraction、LLM judgment。
  - v2 已加入早期 `SRA` 与 `BioProject` provider，用于 raw sequencing context 和项目级上下文补充。
  - 提供 `references/` 中的查询扩展、澄清策略、评分参考和 judgment prompt。

## 如何使用

1. 每个 Skill 目录下都有一个 `SKILL.md`，这是给 AI 助手的运行说明书。
2. `scripts/` 提供稳定脚本能力，`references/` 提供按需读取的规则和参考材料。
3. 与 AI 协作时，优先让助手读取对应 Skill 的 `SKILL.md`，再按需要调用脚本和参考文件。

## 作者

- **Email**: gs.wyuan24@gzu.edu.cn
