# PR #1 独立审查终版（已合并）

## 结论摘要

**最终代码裁定：MERGE；PR #1 已以 merge commit `dd2c358e3458105f31239fe70249445049f47920` 合入 `main`。**

独立核验确认：远端修复分支为 `a0272dd41d19c4a08a45321c5d9adee5627566c9`，远端 `main` 为 `dd2c358...`；该 merge commit 的两个父提交正是旧 main `384b131...` 与 PR head `a0272dd...`；PR head 是 `origin/main` 的祖先，且 `git diff origin/main a0272dd...` 为空，合并没有引入解析差异。

第六轮 P2（`amino_acid:null` 被误分类为 start）已关闭。本地精确 PR head 完整依赖测试 `268/268 OK`、0 skipped；无 `jsonschema` 环境 266 项执行通过、2 项明确 SKIPPED。由于当前到 GitHub Actions API 的 TLS 连接持续 EOF，run `36228676024` 的 success/headSha/268 tests 数据本次无法再次独立读取；这些 CI 数据记录为提交方提供，与已独立确认的远端 refs 和本地测试相容，但不冒充本次独立 API 核验。

合并后仍保留两个非阻塞项：1×P2（start selector 规范化为空仍可绑定无身份 CDS）和 1×P3（null gene/codon 的错误文字误述记录分类规则）。两者均不清除内部 stop ERROR，也不产生自动 PASS；建议合并为一个小型后续 PR。

## 1. 审查范围与版本

- 仓库：`/home/dell/workspace/yuanwei/Bio_project/Bioskills`
- 模块：`mito-assembly-triage/`
- 原 BASE / 旧 main：`384b1317511a510378db1836539134b759fdcf06`
- PR head：`a0272dd41d19c4a08a45321c5d9adee5627566c9`
- merge commit / 当前远端 main：`dd2c358e3458105f31239fe70249445049f47920`
- merge parents：`384b1317511a510378db1836539134b759fdcf06`、`a0272dd41d19c4a08a45321c5d9adee5627566c9`
- `git merge-base --is-ancestor a0272dd origin/main`：exit 0
- `git diff --stat origin/main a0272dd`：空
- BASE→PR head：19 个文件，4184 insertions、189 deletions
- `git diff --check BASE...PR_HEAD`：exit 0
- 离线审查包：`references-review/PR1-FULL-DIFF.patch`，4910 行，SHA-256 `3ddcb3e3bc63383d2114058def200385df77b61991678246cf3e3d225294d4f7`；与 `git diff origin/main-before-merge...a0272dd` 对应内容逐字节一致（合并前已核验）

检查过的全部变更文件：

1. `.github/workflows/community-knowledge.yml`
2. `mito-assembly-triage/SKILL.md`
3. `mito-assembly-triage/references/annotation_quality.md`
4. `mito-assembly-triage/references/diagnostic-playbook.md`
5. `mito-assembly-triage/references/evidence-standard.md`
6. `mito-assembly-triage/references/standard_gene_order.md`
7. `mito-assembly-triage/references/tool-catalog.md`
8. `mito-assembly-triage/scripts/annot_check.py`
9. `mito-assembly-triage/scripts/cox1_id.py`
10. `mito-assembly-triage/tests/fixtures/blast_xml2_outfmt16.xml`
11. `mito-assembly-triage/tests/fixtures/blast_xml_outfmt5.xml`
12. `mito-assembly-triage/tests/gb_fixtures.py`
13. `mito-assembly-triage/tests/test_annot_check_fixes.py`
14. `mito-assembly-triage/tests/test_annot_cox1_edges.py`
15. `mito-assembly-triage/tests/test_annot_semantics.py`
16. `mito-assembly-triage/tests/test_evidence_contracts.py`
17. `mito-assembly-triage/tests/test_pr1_review_regressions.py`
18. `mito-assembly-triage/tests/test_security.py`
19. `mito-assembly-triage/tools/experience.py`

## 2. 审查方法

### 2.1 静态检查

逐轮执行并核对 `git status --short`、HEAD/父提交/BASE、merge-base、完整 diff、逐文件 diff、diff stat、diff check、审查包 hash/cmp、远端 refs、merge parents、祖先关系及合并树差异。

重点审查了：

- partial CDS 的位置对象优先级、正负链、compound/join、跨原点和 5′/3′ 语义；
- `transl_except` 的完整语法消费、位置集合、读框、链方向、aa token、逐 stop 扣除、registry 的 taxon/table/site/selector 绑定；
- BLAST HSP 的 query/subject 重叠、共线性、链方向、coverage/identity 去重、blocked hit、XML2/旧 XML 和故障退出码；
- Schema 的 jsonschema/fallback 等价性、缺失/null、嵌套类型、required/enum/minItems/minLength 和公共 CLI；
- 文档、退出码、科学证据边界与 CI 实际依赖。

`requesting-code-review` 工作流用于固定版本、范围和验收证据；最终结论来自源码、测试、最小复现和 Git 状态，而非开发报告本身。

### 2.2 本地测试

精确 PR head `a0272dd`：

```text
python -m py_compile scripts/*.py tools/*.py tests/*.py
exit 0

python -m unittest discover -s tests -v
Ran 268 tests in 38.592s
OK
```

完整依赖执行 268/268，**0 skipped**。既有 `ResourceWarning` 不影响退出码。

强制无 `jsonschema` 环境：

```text
Ran 268 tests in 34.048s
OK (skipped=2)
```

准确含义：266 项实际执行通过，以下 2 项 **SKIPPED**：

- `SchemaValidatorEquivalenceTests.test_both_paths_agree`
- `SchemaTypeMatrixTests.test_matrix_agrees_with_jsonschema`

合并树与 PR head 内容相同，因此上述独立本地结果适用于合并后的代码树。提交方另报告在 detached `origin/main` worktree 得到相同结果；本报告不把该自报运行替代为独立执行。

### 2.3 CI 核验状态

提交方报告：run `36228676024`，conclusion `success`，headSha=`a0272dd...`，268 tests、0 skipped、diff quality gate 成功。

本次独立复核实际尝试 `gh run view`、GitHub 公共 REST API 和公开页面；CLI/API 均因 TLS `EOF`/`unexpected eof while reading` 无法取回 run 内容。因此精确记录：

```text
远端 refs / merge 状态：独立核验成功
CI run 36228676024 内容：提交方提供；本次独立网络复核 UNAVAILABLE
```

不得把先前 run `36088458880`（针对 `f6e159c`）当作 `a0272dd` 的 CI。

### 2.4 独立最小复现

未修改项目已有测试；复现文件只在 `/tmp`：

- `amino_acid:null`：当前 loader 受控 exit 1，确认第六轮 P2 关闭；
- `amino_acid` 缺键：合法 start 记录仍正常加载；
- `gene:null` / `codon:null`：受控 exit 1，但错误理由存在 P3-1；
- `gene:"?"` start registry + 无身份 CDS + `--tolerate-start :CGA`：实际输出“按已审计例外记录接受”，确认残留 P2-1；
- transl_except 的空 gene、错误 taxon/table/site、无身份 CDS均继续 fail closed，内部 stop ERROR 保留。

## 3. 缺陷清单

### P0

无。

### P1

无。

### P2（已知、非阻塞后续项）

#### P2-1：start registry 和 CLI 仍允许规范化为空的 gene selector

- 文件及准确行号：`scripts/annot_check.py:1136-1139`、`1248-1253`、`803-812`；文档冲突位于 `references/annotation_quality.md:182-186`。
- 触发条件：CDS 没有可识别的 `/gene`/`/product`；registry 的 start 记录 `gene` 为 `"?"`、空字符串或纯标点；CLI 使用 `--tolerate-start :CGA` 或等价空 canonical gene。
- 为什么是真实缺陷：transl_except 分支已经拒绝空 canonical gene，但 start 分支直接写入 `(gene, codon)`，CLI 也只检查 codon 非空。两个空 key 因而能够匹配，并把没有基因身份的记录描述为“按已审计例外记录接受”。文档则笼统声称 selector 归一化为空会加载失败。
- 最小复现：9 bp 无 gene/product CDS `CGA AAA TAA`；registry `{"gene":"?","codon":"CGA",...}`；命令加入 `--tolerate-start :CGA --taxon Lepidoptera`。
- 实际结果：

  ```text
  [WARN ] NONCANONICAL_START_REVIEW: ? 起始密码子 CGA 按已审计例外记录接受;
          taxon=Lepidoptera source=DOI:test rationale=test empty selector
  ```

- 预期正确行为：start registry 和 CLI 均拒绝规范化为空的 gene/codon selector；无身份 CDS 不得绑定 gene-keyed start 证据。
- 最小修改建议：在 start 分支写入前要求 canonical gene 非空、codon 满足三个 IUPAC 碱基；在 `--tolerate-start` 解析处应用相同校验。增加 registry/CLI/无身份 CDS 三层回归测试。
- 严重度依据：会错误描述审计绑定，但非典型起始路径始终是 REVIEW；不会清除内部 stop ERROR，也不会自动 PASS，因此非阻塞。

### P3（非阻塞）

#### P3-1：null gene/codon 的错误信息误述记录分类依据

- 文件及准确行号：`scripts/annot_check.py:1087-1091`
- 触发条件：不含 `amino_acid` 的 start 记录显式设置 `gene:null` 或 `codon:null`。
- 为什么是真实缺陷：拒绝和 exit 1 正确，但消息称“键存在即表示这是 transl_except 记录”；实际上只有 `amino_acid` 键决定类型。
- 最小复现及结果：`{"gene":null,"codon":"CGA",...}` → exit 1，同时输出上述错误理由。
- 预期正确行为：通用消息只说明 selector 不得为 null；仅 `amino_acid` 分支说明键存在决定 transl_except 类型。
- 最小修改建议：拆分错误信息，或删除 gene/codon 情形中的分类说明。

## 4. 已关闭的关键问题

- partial CDS 已以 Biopython fuzzy position 对象为主，冲突显式报告；
- transl_except 已精确绑定读框、链、密码子位置、样本 taxon、遗传密码表和审计 site；未验证 stop 保持 ERROR；
- `amino_acid:null` 不再降格成 start；
- BLAST 缺字段、重叠、非共线、混链和跨原点候选不会进入隐式自动择优；
- query coverage 与 identity 不重复计数；
- Schema fallback 对关键字、类型、显式 null 和嵌套对象的行为与 jsonschema 当前测试矩阵一致；
- 文档不再把注释一致性、BLAST identity 或工程阈值宣称为 reads 支持、物种鉴定或普遍生物学规律。

## 5. 现有测试与 CI 无法覆盖的风险

1. 本次网络故障使 run `36228676024` 无法由审查方再次读取；应保留 run 链接/日志作为项目侧记录。
2. 实时 NCBI QBlast/XML2 未做在线端到端上传验证；现有证据来自协议、本机 Biopython、fixture 与 mock。
3. `source`/`rationale` 的真实性及 `scope="gene_wide"` 是否获文献支持只能人工审计。
4. CI 始终安装 jsonschema；无依赖全套测试需单独环境运行，其中两项交叉比较必然 skipped。
5. Schema 双实现未来仍可能漂移；BLAST 坐标与声明 query/subject 长度的完整边界校验仍是技术债。
6. 新增 null 测试只锁定退出码，没有锁定诊断文字，因而未捕获 P3-1。

## 6. 阻塞与后续处理

PR #1 已合并，且没有遗留 P0/P1。P2-1 与 P3-1 均为非阻塞，但位于同一函数、修复范围很小且测试边界明确。

**后续决定：批准创建一个独立小 PR，同时修复 P2-1 与 P3-1。** 要求：

- 从最新 `main` 建新分支，不复用已合并分支；
- 只修改 registry/CLI selector 校验、诊断文字、对应文档和测试；
- 不借机重构 transl_except、BLAST 或 Schema；
- 覆盖 start registry 与 `--tolerate-start` 的缺键/null/空字符串/纯标点/有效值矩阵；
- 完整测试和无 jsonschema 测试照常运行；
- 新 PR 在 CI 绿灯和独立复核前不合并。

## 7. 最终裁定

**MERGE（已完成，合并内容与获批 PR head 一致）。**

保留意见仅限上述非阻塞 P2/P3 和外部未实测风险，不推翻 PR #1 的合并结论。

本报告归档于 `references-review/PR1-INDEPENDENT-REVIEW.md`。`references-review/` 已被 `.gitignore` 排除，因此这是本地审查档案，不会自动进入后续 PR；除非用户另行明确要求，不应强制加入版本控制。
