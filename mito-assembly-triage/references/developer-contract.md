# 实现契约（格式 / 解析 / 测试约束）

> **本文件是工程约束，不是生物学判据。** 生物证据规则见 `evidence-standard.md` / `annotation_quality.md`；
> 这里的内容回答"代码必须怎样表现、由哪组测试锁住"。
> 拆出本文件的原因：把"selector 归一化后为空则加载失败"这类**实现约束**写在生物学文件里，
> 容易被后来者误读成"领域规则"。

## 1. 案例校验只有一条路径（缺依赖即失败）

案例校验**唯一**的实现是 `schemas/case.schema.json` + `jsonschema`：

| 测试 | 锁住什么 |
|---|---|
| `tests/test_evidence_contracts.py::CaseSchemaVerdictTests` | 固定数据集（合法/缺字段/类型错/非法枚举/嵌套对象错/跨字段冲突，含历史上的分歧用例）的裁定不得改变 |
| `tests/test_pr1_review_regressions.py::SchemaKeywordCoverageTests` | schema 的每个顶层关键字至少有一条“必须被拒绝”的样本 |
| `tests/test_pr1_review_regressions.py::SchemaTypeMatrixTests` | 任意合法 JSON 值都不得让校验器抛异常（只给裁定） |
| `tests/test_missing_dependency.py` | **缺 `jsonschema` 时必须响亮失败**（退出码 3 + 安装命令），且写入口不得落盘 |

**为什么不再有第二份实现**：早期为“没装库的环境”另写了一份手写等价校验
（`_case_errors_without_jsonschema`，88 行）。两份实现都对同一案例说 `VALID`，却靠人肉同步，
并且**已经偏离过两次**（显式 `decision: null`；`case_id` 的 `minLength: 1` 与三个可选数组）。
这违反本 skill 自己的原则——**同一结论必须对应同一条证据路径**；一份可能更弱的校验器
静默给出同一个结论，比直接失败危险得多。因此那份实现被删除（issue #3 以“删除重复实现”结案）。

**缺依赖时的契约**：`case-validate` / `case-anomaly` / `case-reference` 在缺 `jsonschema`
时以**退出码 3** 失败（环境/依赖故障，与 `cox1_id.py` 的 3 同一含义），并打印
`pip install jsonschema`；**不得**降级、不得写入未校验的案例。依赖声明见 `tool_check.md`
与 `scripts/check_env.sh`。

**跨字段业务规则不在 schema 里表达**（如案例级 `RESOLVED` 与某异常 `UNRESOLVED` 并存）；
若将来要加，必须由 schema 表达（必要时用 `allOf`/`if-then`），**不允许**再引入第二份实现。


#### 业务规则层（schema 之外，独立于格式）

`tools/experience.py` 的 `BUSINESS_RULES` 列出一组**跨字段业务关系**，由
`case_business_errors(case, root, verify_inputs=False)` 检查。它们不是格式，因此不写进 schema——
把"格式合法"与"科学自洽"混为一谈，正是上一轮被删掉的那种隐患。每条规则有稳定错误码：

| 错误码 | 含义 |
|---|---|
| `DECISION_ANOMALY_CONFLICT` | 案例级 `decision.status=RESOLVED`，却存在 `UNRESOLVED` 的异常 |
| `ANOMALY_ID_DUPLICATE` | `anomalies[].id` 重复，逐异常判定不可区分 |
| `NORMAL_CASE_HAS_UNRESOLVED_ANOMALY` | `case_type=normal_validation_case` 却挂着未解决异常 |
| `ANOMALY_EVENT_UNKNOWN` | `anomalies[].evidence_events` 引用了 `events.jsonl` 中不存在的事件 |
| `EVENT_HYPOTHESIS_UNKNOWN` | 事件 `impact` 引用了 `hypotheses[]` 中不存在的假设编号 |
| `INPUT_FILE_MISSING` | `inputs[].path` 在磁盘上不存在（**note 级**） |
| `INPUT_SHA256_MISMATCH` | 只有 `--verify-inputs` 时才比对（大 FASTQ 重哈希代价高），不一致即报 |
| `EVENTS_UNPARSABLE` | `events.jsonl` 存在无法解析为 JSON 对象的行 |

**严重度**：`error` 级（记录自相矛盾，如案例判 RESOLVED 却挂 UNRESOLVED 异常）才会让 `case-validate` 判为 INVALID；`note` 级只提示环境性事实（例如归档案例的输入文件不在本机、或长度字段缺失导致某项无法核对），**不改变判读**。输出里每条都标 `[代码/严重度]`。

输出与退出码约定：`case-validate` 把 `FORMAT` 与 `BUSINESS` **分开打印**，`VALID` 只表示
"格式合法 + 已定义业务规则无矛盾"，**不隐含科学正确**；任一层有问题即退出 1，输出里注明是
`INVALID (format)` 还是 `INVALID (business)`。

**写入口不因业务矛盾阻断**：`case-anomaly` / `case-reference` 只**告警**（`⚠ [CODE] …`）并照常写入。
理由：目前没有"修改 decision"的命令，阻断会把用户卡死；格式错误仍然拒绝写入（那时的记录不可用）。

**刻意不作为规则的一项**：`decision.confidence=high` 却没有 `reads_support`。按
`evidence-standard.md` §2.1，无 reads 并不限制与 reads 无关的注释身份结论，因此这**不是**矛盾，
写成规则会与证据标准冲突。此处记录该决定，避免将来被当成漏检。

## 2. `/transl_except` 与例外 registry 的加载约束（fail-closed）

生物学含义（何时算已验证）见 `annotation_quality.md` §9；**格式非法即受控失败**的清单如下：

1. **selector 归一化后为空或显式为 `null` → 加载失败**：`gene` 为 `""`/`"?"`/纯标点/`"(CUN)"`，
   `codon`/`amino_acid` 为空字符串，或键存在但值为 `null`。
   理由：空 selector 会与"缺 `/gene`+`/product` 的 CDS"（canonical `""`）**精确匹配**，
   等于把已审计例外挂到"没有任何基因身份"的记录上。
2. **记录类型按"键是否存在"判定**（而非真值）：`amino_acid` 键一旦出现就是 `transl_except` 记录，
   `""` 或 `null` 均为格式错误，**不得降格**为 start 记录。
3. **`transl_except` 的 `codon`** 必须是三个 IUPAC 碱基；`amino_acid` 必须是合法例外 token（不得为 `TERM`）。
4. **start selector 不得为空**：`--tolerate-start` 的基因名或 start 记录的 `gene` 归一化后为空
   （`":CGA"`、`gene:"?"`/`""`/纯标点）→ 参数/加载受控失败。
5. **统一密码子规则**（start 记录、`transl_except` 记录、`--tolerate-start` 共用同一函数）：
   先 `.strip().upper()`，再要求**恰好三个 IUPAC 碱基**（`[ACGTURYKMSWBDHVN]{3}`）。
   `""` / `"   "` / `"CG"` / `"XXXX"` / `"C1A"` 在加载/参数阶段受控失败；`"cga"` 与 `" CGA "` 归一化为 `"CGA"`。
6. **start 记录按 `(gene, codon)` 唯一**：重复（含仅 `taxon` 不同）→ **加载失败**，不允许后一条静默覆盖；
   若确需多类群共存，必须先扩展 start registry 结构并定义选择规则。
7. **`transl_except` 记录同位点、不同 `taxon`/`transl_table` 可共存**（record list + 显式选择）；
   仅当 selector 全等（gene/codon/amino_acid/位点/taxon/transl_table）才判重复并拒绝。
8. **位点必须绑定**：记录需含 `codon_index` 或 `pos`，或显式 `scope="gene_wide"`；
   `transl_table`/`taxon`/`source`/`rationale` 必须为正确类型且非空。
9. **未登记/未验证只维持 REVIEW**，不升级为已验证：`EXCEPTION_NOT_REGISTERED` / `EXCEPTION_RECORD_INCOMPLETE` /
   `EXCEPTION_TAXON_MISMATCH` / `EXCEPTION_TAXON_UNVERIFIED`；
   而**格式非法（本节第 1–8 项）导致加载失败并退出 1**，不进入 REVIEW，也不做部分生效。

## 3. 退出码语义（不得误读）

| 工具 | 码 | 含义 | **不是** |
|---|---|---|---|
| `annot_check.py` | `0` / `1` / `2` | 无发现 / 有错误（含参数错误）/ 仅待核查 | `2` **不是**通过；`case-validate` 的通过也不是科学验收 |
| `cox1_id.py` | `0` / `1` / `2` / `3` | 得判读 / `insufficient` 或 `no_match` / 拒绝执行 / 网络或结果格式故障 | `3` **不是**"没找到相似序列" |
| `blast_genes.py` | `0` / `1` | 全部基因唯一命中 / 任一基因无唯一命中或命中歧义（**不写部分结果**） | `1` 不是"基因真的缺失" |
| `mitos2_to_genbank.py` | `0` / 非 0 | 写出 `<out.gb>` / 缺文件、坐标越界、`result.fas` 与坐标不一致、未知 feature 类型 | 非 0 时**不得**留下半成品 |
| `depth_analysis.py` | `0` / `1` / `2` | 正常 / 输入或 BAM 错误 / 存在低覆盖区 | `2` 不是"组装错误" |
| `experience.py case-validate/case-anomaly/case-reference` | `0` / `1` / `3` | 通过 / 记录非法或参数错误 / **缺必需依赖 `jsonschema`**（打印 `pip install jsonschema`，不降级） | `3` 不是“案例有问题”，而是校验器不可用 |
| `circularize.py` | `0` / `2` / `1` | 接受候选 / REVIEW（含待回贴）/ 失败或接缝证据不足 | 任何一个码都不能独自证明物理闭环 |

后台任务"子进程正常退出"只表示命令执行成功；科学验收由证据标准另行判定。

### 3.1 经验与案例命令（确切参数）

`SKILL.md` 只给流程与指向（通用参数以 `--help` 为准）；但**涉及证据链记账**的命令在这里给出确切参数，
避免照文档构造出缺必需参数的命令（`--input` 是记录 `inputs[]` + SHA-256 的唯一 CLI 途径，
而 `inputs` 是 `schemas/case.schema.json` 的必需键）。

```bash
# 建档：登记输入文件与哈希（inputs[] 必需），并给出至少一条候选解释
python3 tools/experience.py case-init work/case-001 --issue internal_stop \
  --observation 'nad5 内部 stop' --input assembly_fasta work/asm.fa \
  --hypothesis 'H1: 边界/读码框错误' --hypothesis 'H2: 碱基错误' \
  [--case-type normal_validation_case|tool_failure_case] [--taxon 'Hemiptera: Delphacidae']

# 事件：命令、结果、对假设的影响（如 H1:against）
python3 tools/experience.py case-event work/case-001 --action annot_check \
  --result 'table 5 下仍有内部 stop' --impact H1:against

# 逐异常判定：写入既有 case.json 的 anomalies[]；--event-action 必须命中真实事件
python3 tools/experience.py case-anomaly work/case-001 --id A1 --claim '…' \
  --status UNRESOLVED --confidence low --reads-support NOT_ASSESSED \
  --event-action annot_check [--update]

# 关联已登记的公共参考：id 未登记 / 用途未声明 → 直接失败
python3 tools/experience.py case-reference work/case-001 \
  --reference-id ref-001 --purpose gene_order_comparison [--registry DIR]

python3 tools/experience.py case-validate work/case-001   # 只校验记录格式，不证明科学结论
python3 tools/experience.py case-report  work/case-001   # 生成 case.md 草稿（须按 conclusion-report.md 复核）

# 经验：candidate → verified 的审核链
python3 tools/experience.py propose-lesson --case work/case-001 --next-test '复核边界' \
  [--lesson-domain annotation] [--supporting-case DIR …]
python3 tools/experience.py review-lesson --lesson-id lesson-001 --status verified \
  --reviewer human --reason '记录可核查的独立证据核验'

# 检索与共享（共享/上传是两个独立动作，各自需要授权/校验）
python3 tools/experience.py search-structured --query 'internal_stop nad5'
python3 tools/experience.py export-contribution --case work/case-001 \
  --output contribution.json --authorize
python3 tools/experience.py sync-public --manifest <manifest-url-or-file>
```

代码里的 `required=True` 参数（缺了直接报错）：`case-anomaly --id/--claim/--status/--confidence`、
`case-reference --reference-id/--purpose`、`propose-lesson --case/--next-test`、
`review-lesson --lesson-id/--status/--reviewer/--reason`、`search-structured --query`、
`sync-public --manifest`；`case-init` 必须至少一条 `--hypothesis`。
`--authorize`（导出）与 `--allow-public-upload`（COX1 查询）是**授权闸门**，不是可选装饰。

## 4. 写入安全（工程要求，与生物学无关）

- 候选修复写入**独立目录**，绝不覆盖原始输入；
- 结构化案例的写入是**原子**的（临时文件 + `os.replace`），失败不留 `.tmp`、不留半成品 `case.json`；
- `case-anomaly` 的重复 `id` **fail closed**（文件字节不变，`--update` 才替换）；
- `--event-action` 必须命中 `events.jsonl` 中真实存在的事件（避免悬空关联）；
- `mitos2_to_genbank.py` 的拓扑参数只写声明，**MITOS2 的 circular 模式不是物理环化证据**。

## 4b. 依赖与环境的单一来源

依赖清单只有一个来源：`config/dependencies.json`（tier + probe + purpose + install + stage→requires）。
`tools/env_check.py` 读它做三件事：`--setup` 全量盘点并在 essential 齐全时写
`$MITO_KNOWLEDGE_DIR/environment.lock.json`、`--daily` 轻量检查、`--stage NAME` 单阶段门禁。
`scripts/check_env.sh` 把这三个模式的参数**直通**给该工具，不在 shell 里再抄一份工具清单。

**退出码**（不要混用）：`scripts/check_env.sh` / `env_check.py --setup|--daily` 用 `2` 表示
essential 缺失；`env_check.py --stage` 与**工作脚本**用 `3` 表示"本次所需依赖缺失、该步骤未执行"
（与 `cox1_id.py` 的 3、以及 `case-validate` 缺 `jsonschema` 的 3 同义）。

**脚本门禁**：每个脚本干活前调用 `scripts/_deps.py` 的 `require_stage('<stage>')`，缺依赖即退出 3
（与 `--stage` 同义），**不允许**任何 `except ImportError` 式的降级路径；探测逻辑与 stage→requires
只有一份（`tools/env_check.py` + 清单）。`tests/test_dependency_gates.py` 用"缺 Biopython / 缺 samtools"
两种合成环境钉住它，并断言脚本里不再出现 `except ImportError`。

**不变式**：lock 是缓存，不是信任凭证。日常模式仍真实探测所需工具；只信 lock 就会退化成
"第一次通过、以后永远相信"，即本 skill 反复清除的静默降级。`tests/test_env_check.py` 用
"写 lock 后删掉一个工具，日常模式必须报错"钉住这条。

## 5. 测试与规则的对应关系

失败复现测试**永久保留**；每条硬规则都要有对应的锁：

| 规则 | 锁定测试 |
|---|---|
| GenBank 语义（链、跨原点、`transl_except`、partial） | `tests/test_annot_semantics.py`、`tests/test_annot_check_fixes.py` |
| COX1 结构化 BLAST（XML2/旧 XML/半可读 HSP=格式故障） | `tests/test_annot_cox1_edges.py`、`tests/test_pr1_review_regressions.py` |
| 两条校验路径等价 | `tests/test_evidence_contracts.py`、`tests/test_pr1_review_regressions.py` |
| 工具链桥（MITOS2→GenBank、`run_mitos2` outdir、`check_env` 提示） | `tests/test_toolchain_fixes.py` |
| 案例 CLI（`case-anomaly`、假设编号） | `tests/test_case_anomaly_cli.py` |
| 公共参考登记（accession.version、等级→用途、漂移、截断） | `tests/test_reference_registry.py` |
| 代码 ↔ 文档一致性（flag/错误码/枚举值不得丢文档；规则表可检索） | `tools/audit_code_docs.py` + `tests/test_doc_contracts.py`（基线 `tools/doc_contract_baseline.json`） |
| 案例类型与 lesson 范围/领域（复审核 R-1…R-8） | `tests/test_case_type_and_lesson_scope.py`、`tests/test_review_round_lesson_domain_and_report.py` |
| 注释策略与类群例外 | `tests/test_annotation_policy.py` |
| 数据与路径保护（未发表数据不进公共库） | `.gitignore` 规则 + `git check-ignore`/`git add -n` 人工核验（见 `REAL-DATA` 审计记录） |

新增/修改行为时的顺序：**先写失败复现 → 最小修复 → 补反例 → 跑完整测试与 CI**，CI 绿灯是必要条件。

### 5.1 代码 ↔ 文档一致性审计（常驻护栏）

`tools/audit_code_docs.py` 检查两侧是否漂移，三项任一失败即以非 0 退出（CI 会跑）：

1. **可执行 token 覆盖**：从代码抽取 CLI flag、判定/错误码、枚举值，要求每个 token
   **要么有文档，要么在 `tools/doc_contract_baseline.json` 里显式豁免**。
   新增 token 未写文档也未豁免 → 失败（强制"写文档，或明确豁免"这一次决定）；
   `documented` 列表里的 token 消失 → 失败。
2. **规则表**：一组精选的"代码强制行为"必须仍能在 `SKILL.md` / `references/` 中检索到；
   每条规则自带**代码凭据**，凭据不存在也失败（防止规则表自身腐化）。
3. **基线同步**：`--check-baseline` 要求当前状态与基线完全一致。

**为什么要有它**：一次文档压缩把 6 个 CLI flag 的唯一文档位置删掉了，其中 4 个是
`required=True`；代码没变，所以单元测试全绿——漂移发生在文档侧，需要单独的检查。

**豁免的边界**：可以豁免"逐工具调参 flag"（`--min-mapq`、`--evalue`、`--junction-region` …），
因为 `SKILL.md` 已声明"参数以 `--help` 为准"；但**会改变结论语义或授权语义**的 flag
（`--input`、`--next-test`、`--reference-id`、`--event-action`、`--authorize`、
`--allow-public-upload` 等）必须保持有文档，且已在 `documented` 列表中受保护。

更新基线（`--update-baseline`）是**需要理由的动作**：它意味着"这个 token 故意不写进文档"，
请在 PR 说明里给出理由，而不是用它让 CI 变绿。

## 6. 文档分层约定（避免把实现约束写成领域规则）

| 层 | 归属 | 例子 |
|---|---|---|
| (a) 一般生物学预期 | `annotation_quality.md` / `standard_gene_order.md`（须附来源与覆盖类群） | 13 CDS / 22 tRNA 是**预期不是定律** |
| (b) NCBI 提交审查要求 | `evidence-standard.md` §6（须引用官方页面） | 提交时需向策展人说明差异 |
| (c) 本工具工程阈值 | 相关文件显式标注 | `>8bp` 重叠、tRNA `60–75bp`、`9+/4−`、`max_span_ratio=3.0` |
| (d) 实现约束 | **本文件** | selector 归一化、原子写入、退出码、schema 等价性 |

(c)/(d) 层**不得**写成"NCBI 要求…"或领域公理。

## 7. 已知技术债

1. ~~**两套案例校验实现**~~ **已解决**：手写兜底实现已删除，改为单一 `jsonschema` 路径 +
   缺依赖时以退出码 3 失败（见 §1）。固定数据集测试保留，继续钉住裁定。

2. **跨字段矛盾不校验**（如案例级 `RESOLVED` 与异常 `UNRESOLVED` 并存）：schema 不表达，
   两套实现都不拒绝；若要加必须同时加在两边。
3. **BLAST 坐标边界**：`cox1_id.py` 的 HSP 边界情形仍有未覆盖组合（详见 issue 记录）。
4. **`case_type` 是可选字段**：`learning-policy.md` 把记录分类列为必需，但 schema 未列入 `required`。
   若改成强制项，已有/外部归档案例（可能缺该字段）会一次性变成 `INVALID`，
   因此当前策略是**可选 + 缺失时按 `abnormal_case` 保守处理 + 新案例一律写入**；
   要真正强制，必须同时提供迁移方案（回填 `case_type`）与两套校验路径的同步修改。
5. **lesson 的领域/范围不是 schema 强制项**：写入端（`propose-lesson`）是 fail-closed 的，
   验收端（`validate_public_lesson`）对缺失值做保守补齐后才接受；两边规则不同步就会再造出偏差，
   因此改任意一边必须同时跑 `tests/test_review_round_lesson_domain_and_report.py`。

6. **报告按状态机械分类**：`case_report()` 将 RESOLVED/NO_CHANGE 放入“有证据支持”，
   UNRESOLVED 放入“无证据支持”。这不符合证据语义；当前只能把输出当草稿，按
   [conclusion-report.md](conclusion-report.md) 复核。下一步修复应以部分支持但未解决、
   证据冲突、无证据的 NO_CHANGE 等案例验证，不能用改标题掩盖证据欠缺。
7. **科学报告字段未完整承载**：schema/CLI 未完整校验支持/反证、未测项、范围、局限、
   结论类型与复核信息。`--event-action` 是可选的，只验证被传入的 action 存在；
   action 同名时无法唯一定位某次执行。需由事件与报告补足，不能称“格式通过即完整可追溯”。
   `case-anomaly --update` 会替换该条记录；手工附加字段不能依赖它保留。
8. **类群配置与经验独立性**：`--taxon` 不自动切换昆虫阈值；lesson 的“不同 case_id
   且不同类群”是当前推广门槛的工程计数，不等于科学独立性。文档解释层不得冒称代码已修复。
9. **覆盖提示过强**：`depth_analysis.py` 会将高 soft-clip 比例按区域标签提示为“正常长度异质性”
   或“缺失序列”。这些是未验证的解释，不能作为诊断事实；当前按证据标准重新判读。

## 8. 注释解析细节（仅调试或维护时读取）

以下保留原注释手册中的解析/格式约束；科学接受条件仍见
[annotation_quality.md](annotation_quality.md) §2/§9。下面的 JSON 仅演示格式，含占位来源，
不是可直接使用的已审计生物学例外。

   - `/transl_except` **分三层，存在不等于豁免**：
     **① 语法层**：代码按 CDS 转录顺序建立每个密码子的**精确基因组位置集合**；整个 qualifier 必须被
     **完整消费**（一个或多个 `(pos:...,aa:...)`），语法不完整就记 `TRANSL_EXCEPT_UNPARSED` 并**整条作废，
     不沿用合法前缀**（尾随逗号、尾随损坏文本、缺 `pos:` token、模糊位置、括号不配对都属此类）。
     `pos` 支持 `a..b`、`complement(a..b)`、`join(..)`/`order(..)`（密码子可跨 `join()` 边界）；
     同一 `pos` 内**不同链方向的片段混用直接拒绝**（整段位置不得由 OR 合并成单链）。
     `aa` 须是合法三字母代码且**不是 `TERM`**；`pos` 方向须与 CDS 一致（**负链 CDS 必须 `pos:complement(a..b)`**，
     正链不得写 `complement`）；位置集合必须**恰好**等于某个真实内部终止密码子。通过只记
     `TRANSL_EXCEPT_MATCHED`（位置/读框事实），**MATCHED 本身不等于已接受**。
     **② 生物学层**：语法匹配后还必须有审计证据才能接受。证据来自 `--exception-registry` 的
     `transl_except` 条目，键为 `gene + codon + amino_acid`，且必须满足三个 fail-closed 约束：
     （a）**样本绑定**：必须显式给出 `--taxon` 并与记录 `taxon` 一致 —— 未提供 `--taxon` 时无法把记录绑定到本样本，
     **不得**升为已验证；
     （b）**位点绑定**：记录必须说明它覆盖**哪个** stop —— 给 `codon_index`（CDS 转录顺序上的 1-based 密码子序号）
     或 `pos`（该密码子的基因组位置，如 `10..12`）之一；只有**显式** `scope="gene_wide"` 的记录才允许复用多个位点。
     不绑定位点的记录在**加载时**受控失败（否则一条记录会静默变成基因范围的重编码规则）；
     （c）**字段完整且类型正确**：`transl_table` 须为正整数且等于 `--table`；`source`、`rationale` 须为非空字符串；
     数组/对象/整数等错误 JSON 类型在**加载时**受控失败（`str(value)` 非空不等于已审计）。
     三者均满足才记 `TRANSL_EXCEPT_VALIDATED`（并输出 `scope=`）并扣除该项。
     **不引入全局密码子→氨基酸重编码表**：密码子含义随类群与密码表变化，把“合法 INSDC token”当成
     “已验证例外”会把类群特异的例外伪装成普遍规律。因此 `TAA -> Gln` 这类任意声明即使位置匹配，
     也只能得到 `TRANSL_EXCEPT_DECLARED_UNVERIFIED`（REVIEW），**内部终止仍保留 ERROR**。
     **③ 结果层**：因此**一条声明不能吸收两个内部终止**；一条未绑定该位点的记录也不能验证别处的同类 stop；
     未验证、未解释的内部终止仍是 ERROR。
     registry 的最小 transl_except 记录：
     `{"gene":"cox1","codon":"TAA","amino_acid":"Trp","pos":"10..12",`
     `"transl_table":5,"taxon":"Lepidoptera","source":"DOI ...","rationale":"..."}`
     （或用 `"codon_index":4`，或 `"scope":"gene_wide"` 且不给位点）。

partial 的实现约束：

   - 代码直接读 Biopython 的 `BeforePosition`/`AfterPosition` **位置对象**（不解析 location 字符串），
     按链方向解释生物学的 5'/3' 端；**负链的 5' 端在高坐标**，跨原点 `join()` 的各段按转录顺序排列；
     位置对象**始终优先于**字符串兜底（仅当 parts 完全不带 fuzzy 信息时才用字符串），
     两者矛盾时以对象为准**并报** `PARTIAL_SOURCE_CONFLICT`，不静默取值；

## 9. COX1 BLAST 结果解析契约

查询片段为 400–5000 bp，结果只提供分类线索，不能单独确定物种。它请求 `FORMAT_TYPE=XML2`，并**同时**能读
XML2（`-outfmt 16`）与旧版 XML（`-outfmt 5`）两种方言；`FORMAT_OBJECT=SearchInfo` 按 NCBI 实际返回的
QBlast 文本（`RID =` / `Status=`）解析，XML 形式也兼容。**每个 `<Hsp>` 的 query/hit 坐标、identity、align-len、
bit-score 都必须存在且数值合法**：缺任一必需字段时整份结果按**格式故障**（退出码 3）处理 ——
缺少 subject 坐标就无法做方向/共线性/目标重用检查，不能靠“只剩单个 HSP”绕过验证。它报告每个候选的 **query coverage**（多 HSP 并集）与 identity，并把 RID 轮询与结果下载分开。
退出码把任务状态与判读状态分开：`0` 得判读 / `1` insufficient 或 no_match（分析结论）/ `2` 拒绝执行 /
`3` 网络或结果格式故障。多 HSP 的指标只有同时满足**不重叠**与**目标共线**才回总，检查分开做：
① 全部 HSP 必须同一目标链（真正的正负链混合属矛盾证据）；
② 按 `query_from` 排序后，`hit_from` 必须**沿目标链方向单调推进**（正链递增、**负链递减**）；
③ **query 侧与 target 侧都不能复用同一批碱基**：query 区间重叠记 `overlapping_hsps`，target 区间重叠记
`subject_overlap`（两条 query 分别打到目标同一段的重复/塌缩情形，仍不是两条独立比对）；
其中 **`max_span_ratio=3.0`（目标跨度 ≤ 3 倍比对长度）是额外的工程预警，不是 COX1 生物学标准，也不能代替方向检查**。
坐标约定经本地真实 `blastn` 实测锁定：连续负链 `q1..400→s1800..1401`、`q501..900→s1400..1001`（应**接受**）；
同两段但顺序倒置 `q1..400→s1400..1001`、`q500..900→s1801..1401`（应**拒绝**）。
回总失败时记 `overlapping_hsps` / `subject_overlap` / `non_collinear_hsps` + `CONFLICTING_ALIGNMENT`；若跳变同时触及目标两端，
另记 `CROSS_ORIGIN_CANDIDATE` —— 环状参考下可能是跨原点排列，但**必须**有明确坐标与结构证据，不得直接按连续线性比对接受。
这些情形都**不参与自动择优**，但 `AMBIGUOUS_ALIGNMENT` 的含义是“**多 HSP 指标无法可靠汇总**”，
**不是**“该 hit 不是有效候选”：CLI 会逐条打印被阻断候选的 HSP 明细（`q.. → s..  identity bits`），
并用 `--output-json` 把每个候选的全部 HSP 与 blocker 持久化；该文件无论判读结果如何都会写出。最优 accession 不等于已完成物种鉴定。
