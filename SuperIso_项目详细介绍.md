# SuperIso 项目详细介绍

## 一、项目定位与物理目标

SuperIso 是一个面向味物理与精密观测量分析的数值计算框架，其核心用途不是直接生成新物理模型的质量谱，而是在给定标准模型或新物理模型参数点之后，系统计算一大批低能观测量，并进一步判断这些理论预测与实验数据之间是否一致。从这个意义上说，SuperIso 处在“模型参数”和“实验限制”之间的中间层：上游可以是标准模型、MSSM、NMSSM、2HDM 或者某个已经写成 SLHA 文件的有效参数点，下游则是各种 flavour observables、稀有衰变分支比、角分布观测量、轻子味道普适性比值以及全局 $\chi^2$ 统计量。项目总览写在 [README.md](README.md) 中，其中明确说明该程序用于计算 flavour physics observables、muon anomalous magnetic moment，并支持通过 Isajet、Softsusy、SPheno、SuSpect、NMSSMTools 与 2HDMC 等外部工具得到的模型参数点。

从物理研究的角度看，这类程序的价值在于，许多稀有衰变过程在标准模型中本来就是经过回路抑制或 CKM 抑制的，因此理论预测往往较小，而新物理若在相同有效算符上产生附加贡献，即使修正幅度不大，也可能对观测量造成可见偏移。也正因为如此，味物理观测量经常能对新物理模型给出比直接搜索更强、或者至少非常互补的间接限制。SuperIso 的目标正是把这种“从模型到观测量，再到限制”的链条自动化和模块化。

## 二、项目整体结构与代码组织

这个项目的代码结构可以分成三层理解。第一层是根目录中的主程序文件，它们定义了用户从什么输入出发、以什么形式得到输出。例如，程序 [sm.c](sm.c) 用于在标准模型下直接计算若干典型观测量，[slha.c](slha.c) 用于读取用户给定的 SLHA 文件并对该参数点进行 flavour 计算，[sm_chi2.c](sm_chi2.c) 和 [slha_chi2.c](slha_chi2.c) 则在预测值之外进一步加入实验输入与协方差矩阵，从而给出全局 $\chi^2$。除此之外，根目录中还有一组专门针对不同新物理场景的入口，例如 [cmssm.c](cmssm.c)、[amsb.c](amsb.c)、[gmsb.c](gmsb.c)、[nuhm.c](nuhm.c)、[thdm.c](thdm.c)、[cnmssm.c](cnmssm.c) 等，这些程序并不重新实现 flavour 物理，而是负责为不同模型构造参数点，再调用底层通用计算模块。

第二层是 `src/` 目录中的核心物理计算模块，这一层才是真正实现 flavour observable 的地方。按照功能划分，这些源文件大体可以分成参数与输入模块、Wilson 系数与 RGE 模块、具体观测量模块，以及统计分析模块。比如 [src/leshouches.c](src/leshouches.c) 负责参数结构初始化、SLHA 相关设置和若干理论选项；[src/wilson.c](src/wilson.c) 负责 Wilson 系数与其演化；[src/bsgamma.c](src/bsgamma.c)、[src/bkstarll.c](src/bkstarll.c)、[src/bkll.c](src/bkll.c)、[src/bsll.c](src/bsll.c)、[src/bsmumu.c](src/bsmumu.c)、[src/kaon.c](src/kaon.c) 等分别负责不同衰变道的理论计算；而 [src/chi2.c](src/chi2.c) 则负责不确定度、协方差矩阵和全局拟合。

第三层是 `chi2_input/` 目录中的数据文件。这些文件不是程序逻辑本身，却是全局拟合不可缺少的一部分。像 [chi2_input/myobs.in](chi2_input/myobs.in) 给出了默认参与拟合的观测量列表，[chi2_input/exp_values.in](chi2_input/exp_values.in) 存放实验中心值和误差，[chi2_input/exp_corr.in](chi2_input/exp_corr.in) 给出实验相关矩阵，[chi2_input/nuisance.in](chi2_input/nuisance.in) 和 [chi2_input/nuisance_corr.in](chi2_input/nuisance_corr.in) 则定义理论 nuisance 参数及其相关性。也就是说，SuperIso 不仅是一套公式实现，它还包含了一套可直接用于 flavour global fit 的输入数据库。

## 三、程序从输入到输出的完整计算流程

SuperIso 的一个典型执行流程，可以用 [slha_chi2.c](slha_chi2.c) 来概括，因为这个程序几乎经过了项目中最完整的所有步骤。首先，程序读取用户给定的 SLHA 文件名，并通过 `test_slha` 判断该文件是否存在、格式是否正确、模型是否在当前版本中受支持。接着，程序从 [chi2_input/myobs.in](chi2_input/myobs.in) 中读取要比较的观测量名称列表，再从 [chi2_input/exp_values.in](chi2_input/exp_values.in)、[chi2_input/exp_corr.in](chi2_input/exp_corr.in)、[chi2_input/nuisance.in](chi2_input/nuisance.in) 和 [chi2_input/nuisance_corr.in](chi2_input/nuisance_corr.in) 中装入实验数据和理论 nuisance 输入。随后，程序初始化 `parameters` 结构体，把 SLHA 中的参数读入内部数据结构，再基于这一组参数计算所有观测量的理论预测值。

在预测值得到之后，程序并不会立刻把理论值和实验中心值逐项作差，而是继续构造理论协方差矩阵和实验协方差矩阵。理论协方差用于描述 CKM 元素、强耦合常数、夸克质量、衰变常数、形状因子和幂修正等理论不确定度如何共同影响各个观测量；实验协方差则描述实验测量之间的统计和系统相关性。总协方差矩阵由这两部分相加得到，而全局统计量则由标准二次型

$$
\chi^2 = (O_{\mathrm{th}}-O_{\mathrm{exp}})^T V^{-1} (O_{\mathrm{th}}-O_{\mathrm{exp}})
$$

给出，其中 $O_{\mathrm{th}}$ 是理论预测向量，$O_{\mathrm{exp}}$ 是实验中心值向量，$V$ 是总协方差矩阵。这一套流程在 [slha_chi2.c](slha_chi2.c) 和 [sm_chi2.c](sm_chi2.c) 中都能清楚看到，因此可以说，SuperIso 的统计比较是建立在向量与矩阵层面的全局分析，而不是对单个观测量逐项手工比较。

如果用户并不需要做全局拟合，而只想看某一个参数点的代表性观测量，那么 [slha.c](slha.c) 就提供了较简洁的流程。该程序在读入 SLHA 文件后，会先判断该点是否有效，然后依次输出常用的 flavour observable，例如 $\mathrm{BR}(b\to s\gamma)$、$\mathrm{BR}(B_s\to\mu^+\mu^-)$、$B\to K^*\mu^+\mu^-$ 低、高 $q^2$ 区间中的 $\mathrm{AFB}$、$F_L$、$P_i'$ 等量，并生成 FLHA 格式的输出文件 [output.flha](output.flha)。这说明 SuperIso 同时支持“快速看一个点”和“把一组点拿去做全局统计”两种使用场景。

## 四、参数初始化、理论选项与形状因子设定

项目中最基础的参数设置由 [src/leshouches.c](src/leshouches.c) 中的 `Init_param` 完成。这个函数不只是把结构体成员清零，更重要的是它定义了半轻子衰变计算所采用的理论方案。例如，程序默认设置 `param->fullFF=1`，表示在半轻子衰变中使用 full form factor approach；同时还定义了 `BKstar_implementation` 来选择 $B\to K^*\ell\ell$ 中幂修正和 hadronic contribution 的实现方式。对于 $B\to K^*$、$B_s\to\phi$ 和 $B\to K$ 的形状因子，程序还提供了多个不同文献来源的参数化方案，分别通过 `BKstar_FormFactor_choice`、`Bsphi_FormFactor_choice` 和 `BK_FormFactor_choice` 控制。

这一点的物理意义非常重要。因为 flavour observable 的理论预测不只依赖 Wilson 系数，还依赖形状因子、非局域 charm loop、幂修正和低能强子物理的处理方式。不同理论输入方案虽然不改变模型本身，但会改变预测值的中心值和误差带，进而影响参数拟合结果的稳定性。因此，SuperIso 在参数初始化层面就把这些理论自由度显式暴露给用户，而不是把它们硬编码成不可更改的内部常数。

## 五、Wilson 系数、有效哈密顿量与 RGE running

整个项目的 flavour physics 核心建立在有效哈密顿量方法之上。也就是说，无论上游输入是 SM 还是某个新物理模型，程序最终都要把模型信息转化为低能有效算符的 Wilson 系数，然后再由这些 Wilson 系数去决定具体观测量。README 对 [src/wilson.c](src/wilson.c) 的描述就是 “calculation of the Wilson coefficients and RGE runnings”，这已经点出了该文件的核心角色。

在 [sm.c](sm.c) 中可以直接看到这条计算链的典型写法。程序首先定义一组复数数组来保存不同标度下的 Wilson 系数，例如 `C0w`、`C1w`、`C0b`、`C1b`、`C2b` 和 chirality-flipped 系数 `Cpb`。随后，代码依次调用 `CW_calculator`、`C_calculator_base1`、`C_calculator_base2` 和 `Cprime_calculator` 等函数，先在高标度构造系数，再把它们演化到低标度 $\mu_W$、$\mu_b$ 或 spectator 相关的标度上，最后再把这些系数送入具体的观测量模块，例如 `bsgamma`、`delta0`、`BR_BKstargamma` 和 `Bsmumu_untag`。这说明从代码实现上看，程序并不是把每个 flavour observable 当作黑箱，而是严格遵循“高能匹配 $\to$ RGE 演化 $\to$ 低能矩阵元 $\to$ 观测量”的有效理论逻辑。

这种实现方式的优点在于，不同模型之间可以共享同一套 flavour observable 模块。对于 2HDM、MSSM、NMSSM 等模型而言，上游模型依赖主要体现在 Wilson 系数的来源不同；而一旦这些系数在低能标度上被确定，后续的 flavour 计算流程在很大程度上就是统一的。这也是 SuperIso 能同时支持多种模型场景的根本原因。

## 六、项目中包含的观测量类别

如果只看 [sm.c](sm.c) 和 [slha.c](slha.c)，用户会看到一组常用 flavour observable 的打印输出，例如 $\mathrm{BR}(b\to s\gamma)$、$\delta_0(B\to K^*\gamma)$、$\mathrm{BR}(B_s\to\mu^+\mu^-)$、$B\to K^*\mu^+\mu^-$ 在低、高 $q^2$ 区间中的 $\mathrm{AFB}$、$F_L$、$P_1$、$P_2$、$P_4'$、$P_5'$、$P_6'$、$P_8'$ 等，以及 $B\to\tau\nu$、$B\to D^{(*)}\tau\nu$、$D_s\to\ell\nu$、$D\to\mu\nu$、$K\to\mu\nu$、$K\to\pi\nu\nu$、$K_L\to\ell\ell$ 和 muon $g-2$。这些常用输出足以覆盖 flavour phenomenology 中最常见的观测量类别，因此对于扫描参数点和快速判断某个模型是否大体可行已经非常有用。

但是，项目真正支持的观测量范围远不止这些。最完整的列表要从 [src/interpreter.c](src/interpreter.c) 与 [chi2_input/myobs.in](chi2_input/myobs.in) 两个地方一起看。前者定义了支持的观测量“类型”和“衰变道”，后者给出了当前全局拟合默认纳入的具体 observable 名称。根据 [src/interpreter.c](src/interpreter.c)，程序支持的观测量类型包括 $\mathrm{BR}$、$\mathrm{BRuntag}$、$d\Gamma/dq^2$、$R-1$、$\mathrm{AI}$、$\mathrm{ACP}$、$\mathrm{AFB}$、$F_L$、$F_T$、$F_H$、$A_T^{(i)}$、$H_T^{(i)}$、$S_i$、$P_i$、$P_i'$ 以及多种 CP 版本的角观测量，还允许直接把 Wilson 系数、其虚实部、CKM 元素、甚至形状因子本身作为输出对象。这意味着 interpreter 实际上是一个通用观测量调度层，能够把字符串形式的 observable 名称映射成具体函数调用。

从 [chi2_input/myobs.in](chi2_input/myobs.in) 可以看到，默认参与全局拟合的观测量覆盖了若干大类。第一类是辐射和纯轻子稀有衰变，例如 `BR_BXsgamma`、`AI_BKstargamma`、`BRuntag_Bsmumu` 等。第二类是包容和专属的 $b\to s\ell\ell$ 过程，例如 `BR_BXsmumu_1_6`、`BR_BXsee_14.2_22`、`dGamma/dq2_BKstarmumu_1.1_6`、`FL_BKstarmumu_0.1_0.98` 等。第三类是角观测量，例如 `S3_BKstarmumu_1.1_2.5`、`P4prime_B0Kstar0mumu_1.1_2_CMS`、`P5prime_B0Kstar0mumu_2_4.3_CMS` 等，它们反映的是不同 transversity 振幅之间的干涉信息。第四类是轻子味道普适性比值，例如 `R-1_BKll_0.1_1.1`、`R-1_BKstarll_0.045_6` 和 `R-1_Bsphill_1.1_6` 等。第五类是更复杂的专属衰变，如 $B_s\to\phi\mu\mu$ 和 $\Lambda_b\to\Lambda\mu\mu$ 的微分衰变率与角分布观测量。总的来说，SuperIso 并不是只覆盖几项“代表性 observable”，而是提供了一套非常完整的 flavour global fit 基础设施。

## 七、各类观测量在程序中的具体实现方式

### 1. 辐射衰变：$b\to s\gamma$ 与 $B\to K^*\gamma$

在 [src/bsgamma.c](src/bsgamma.c) 的开头可以看到大量以 `phi77`、`phi78`、`phi22`、`phi27`、`phi47` 等命名的函数。这些函数对应的并不是简单的数值插值，而是辐射衰变率计算中来自不同算符、不同微扰阶以及不同相空间部分的核函数和修正项。从结构上看，程序首先通过 Wilson 系数模块得到诸如 $C_7$、$C_8$ 及四夸克算符的低能系数，再把这些系数与 [src/bsgamma.c](src/bsgamma.c) 中定义的核函数组合起来，从而得到包容型 $b\to s\gamma$ 的衰变率。对于 $B\to K^*\gamma$ 和同位旋不对称 $\delta_0(B\to K^*\gamma)$，程序还会额外考虑 spectator 相关标度和强子输入，因此在 [sm.c](sm.c) 中除了 $\mu_b$ 之外还定义了 `mu_spec`，并调用了 `C_calculator_base2` 来在不同标度下评估所需的 Wilson 系数。

物理上，这些观测量的约束力来源于它们对电磁企鹅算符尤其敏感，特别是 $C_7$ 和 $C_8$。任何会显著改变电磁偶极算符系数的新粒子回路，例如 2HDM 中的带电 Higgs、MSSM 中的 chargino-stop 回路或 gluino-squark 回路，都可能在这些观测量中留下明显痕迹。因此，$\mathrm{BR}(B\to X_s\gamma)$ 往往是限制 2HDM 和 SUSY 参数空间最强的 flavour observable 之一。

### 2. 包容型 $b\to s\ell\ell$

在 [src/bsll.c](src/bsll.c) 中，程序实现了包容型 $B\to X_s\ell^+\ell^-$ 的一整套函数。例如 `g_bsll`、`sigma7_bsll`、`sigma9_bsll`、`f7_bsll` 和 `f9_bsll` 等函数，体现的是微扰 QCD 修正、不同 Wilson 系数的组合以及短程贡献的结构。更值得注意的是，代码中还存在 `Rcchad` 和 `g_bsll_parametrized` 这样的函数，用于参数化 charmonium resonance 与连续谱贡献，这说明作者没有把 inclusive $b\to s\ell\ell$ 简化成纯 partonic 近似，而是显式考虑了 charm resonance 对低能区域的影响。

这类观测量的重要性在于，它们对 $C_7$、$C_9$、$C_{10}$ 等有效算符都敏感，并且不同 $q^2$ 区域的敏感性并不相同。低 $q^2$ 区域常常对光子极点和 $C_7$ 更敏感，而高 $q^2$ 区域则更多受 $C_9$、$C_{10}$ 及非局域 charm 贡献影响。因此，程序中把不同 $q^2$ 区间分开处理，并在 [chi2_input/myobs.in](chi2_input/myobs.in) 中用不同 bin 的 observable 名称明确区分，这是非常自然且必要的做法。

### 3. 专属半轻子衰变：$B\to K^*\ell\ell$、$B\to K\ell\ell$、$B_s\to\phi\ell\ell$

这部分是项目中最复杂、也是物理信息最丰富的一块。在 [src/bkstarll.c](src/bkstarll.c) 中，可以看到诸如 `h_bkll`、`phi_Kstar`、`I1_bkll`、`tperp_bkll`、`tpar_bkll`、`A_Seidel`、`B_Seidel` 等函数。仅从命名和结构就可以看出，这些实现围绕着 helicity 或 transversity 振幅展开，并包含了 charm loop、非局域贡献、形状因子卷积以及不同极化分量的处理。最终，这些振幅被重组为实验上真正测量的观测量，比如微分衰变率 $d\Gamma/dq^2$、前后向不对称 $A_{\mathrm{FB}}$、纵向极化分数 $F_L$，以及优化构造的 $P_i$ 和 $P_i'$ 观测量。

这一类观测量之所以在近年来 flavour anomaly 研究中异常重要，是因为它们不仅对总分支比敏感，更对不同算符之间的干涉结构敏感。简单地说，若新物理改变了 $C_9$，那么它对 $P_5'$、$A_{\mathrm{FB}}$、$F_L$ 和微分谱的影响模式通常不同于改变 $C_{10}$、$C_7$ 或右手流算符时的模式。因此，通过比较一整套角观测量而不是单个分支比，可以更准确地区分新物理的算符结构。SuperIso 用分 bin 的方式把这些观测量组织在 [chi2_input/myobs.in](chi2_input/myobs.in) 中，并通过 [src/interpreter.c](src/interpreter.c) 统一调度它们的计算，因此它在代码架构上天然适合做这类全局 EFT 分析。

### 4. 纯轻子衰变：$B_s\to\mu^+\mu^-$ 与相关过程

在 [sm.c](sm.c) 和 [slha.c](slha.c) 中，$\mathrm{BR}(B_s\to\mu^+\mu^-)$、$\mathrm{BR}(B_s\to\mu^+\mu^-)_{\text{untag}}$ 和 $\mathrm{BR}(B_d\to\mu^+\mu^-)$ 都是直接输出的基本观测量。这些过程在理论上相对洁净，因为末态简单，强子输入主要集中在衰变常数上，所以它们非常适合用来检验轴矢量、标量和赝标量型有效算符。特别是在大 $\tan\beta$ 的 SUSY 或扩展 Higgs 场景中，标量交换可能显著增强这类衰变，因此它们对 Higgs sector 的约束往往很强。程序中同时给出 tagged 与 untagged 版本，也是因为实验测量与理论定义之间在时间积分等细节上存在差异，完整处理需要显式区分。

### 5. 带电流衰变：$B\to\tau\nu$、$B\to D^{(*)}\tau\nu$、$D_s\to\ell\nu$、$K\to\mu\nu$

这部分观测量在 [sm.c](sm.c) 中有非常清晰的展示。程序不仅输出总分支比，还输出如 $A_{\mathrm{FB}}$、$P_\tau$、$P_{D^*}$ 以及不同轻子模式之间的比值。与稀有 FCNC 过程不同，这些衰变在标准模型中是树级过程，因此它们不是通过 loop 抑制体现出新物理敏感性，而是通过新标量流、右手流或额外 charged current 结构改变轻子质量依赖性和角分布结构。比如带电 Higgs 就可以直接修改与 $\tau$ 相关的衰变幅，从而让 $R(D)$、$R(D^*)$、$R(B\to\tau\nu)$ 等量偏离标准模型预测。也正因为如此，这类观测量通常和稀有 $b\to s$ 过程形成互补约束：前者更敏感于 charged current 和标量耦合，后者更敏感于 FCNC 和企鹅算符。

### 6. Kaon 观测量与 muon $g-2$

在 [sm.c](sm.c) 的后半部分，程序还计算了 $K^+\to\pi^+\nu\bar\nu$、$K_L\to\pi^0\nu\bar\nu$、$K_L\to\mu^+\mu^-$、$K_L\to e^+e^-$、$K_S\to\mu^+\mu^-$、$K_L\to\pi^0 e^+e^-$ 和 $K_L\to\pi^0\mu^+\mu^-$ 等过程。Kaon 体系在 flavour physics 中历来重要，因为它对 $s\to d$ 的短程物理和 CP 结构极其敏感，而标准模型预测通常较为干净，因此很适合测试额外 flavour violation 和新的 CP phase。至于 muon $g-2$，虽然它严格来说不属于 flavour observable，但由于它同样是一个对新粒子回路极为敏感的精密量，而且和 flavour observable 常常同时被用来约束同一模型参数空间，所以在 SuperIso 中被自然地纳入了同一套分析工具链。

## 八、interpreter 层的作用：为什么观测量可以用字符串统一管理

SuperIso 中一个非常精巧的设计是 [src/interpreter.c](src/interpreter.c)。在大多数 flavour 程序里，不同观测量常常由不同函数分别处理，用户若想添加或选择 observables，往往必须改代码。而在这个项目里，作者通过 interpreter 的方式，把“观测量名称”抽象成统一的输入接口。例如 `BR_BXsgamma`、`FL_BKstarmumu_1.1_2.5` 或 `P5prime_B0Kstar0mumu_2_4.3_CMS` 这类字符串，会先被拆解成 observable type、衰变道和 bin 范围，然后再由 `compute_nameobs` 调用相应模块给出数值。

这一设计的直接好处是，统计程序不需要关心每个 observable 的具体来源。它只需要从 [chi2_input/myobs.in](chi2_input/myobs.in) 逐行读取 observable 名称，交给 interpreter，再由 interpreter 统一返回预测值。这样一来，往拟合中增加一个新的观测量，通常不需要重写整体统计框架，只需要保证该 observable 的字符串名字能够被 interpreter 正确识别，并且对应的底层物理模块已经实现。这种架构使得项目具有很好的可扩展性，也使得 flavour global fit 的输入组织方式更加清晰。

## 九、这些观测量为什么能对模型产生限制

从物理上讲，SuperIso 中这些观测量之所以有约束力，主要有三层原因。第一，许多 flavour process 在标准模型中本身就被 GIM 机制、CKM 结构或 loop 结构强烈抑制，因此即使新物理只带来小幅修正，也可能与标准模型贡献同量级。第二，不同观测量对 Wilson 系数的依赖不同，有的主要敏感于 $C_7$，有的主要敏感于 $C_9$ 和 $C_{10}$，有的则对标量和赝标量算符更加敏感，因此通过同时比较多类 observables，可以较清楚地判断哪一类算符结构可能被新物理改变。第三，某些观测量的理论不确定度相对较小，或者可以通过构造优化观测量和比值来减少对形状因子的依赖，从而增强对新物理信号的分辨能力。

以辐射衰变为例，$\mathrm{BR}(B\to X_s\gamma)$ 对电磁偶极算符非常敏感，因此任何显著改变 $C_7$ 或 $C_8$ 的模型都会受到强限制。以 $B_s\to\mu^+\mu^-$ 为例，它对轴矢量和标量算符特别敏感，所以扩展 Higgs 扇区或大 $\tan\beta$ SUSY 场景往往会受到它的强限制。以 $B\to K^*\mu^+\mu^-$ 的角观测量为例，诸如 $P_5'$、$P_4'$、$A_{\mathrm{FB}}$ 和 $F_L$ 等量，能够帮助区分新物理究竟更像是改变了 $C_9$、$C_{10}$、右手流，还是主要表现为 hadronic uncertainty。至于 $R(K)$、$R(K^*)$ 一类轻子味道普适性比值，它们之所以极具威力，是因为在标准模型中电子和缪子通道本应高度一致，若出现显著偏离，就很自然地指向轻子非普适的新相互作用。

## 十、这些限制主要针对哪些模型

从项目入口程序和 README 的说明可以看出，SuperIso 的 flavour observable 可用于限制多类新物理模型。对于 2HDM 而言，最典型的约束来自 $b\to s\gamma$，因为带电 Higgs 会直接进入电磁企鹅图并改变 $C_7$。这也是为什么 2HDM 的 flavour phenomenology 中，$\mathrm{BR}(B\to X_s\gamma)$ 常常给出对带电 Higgs 质量最强的下限之一。对于 MSSM 及其受限场景，如 CMSSM、AMSB、GMSB 和 NUHM，flavour observable 对 chargino-stop、gluino-squark、Higgs sector 以及 slepton/electroweakino 扇区都可能形成约束。比如 $b\to s\gamma$ 对有色超对称粒子和带电 Higgs 的回路很敏感，$B_s\to\mu^+\mu^-$ 对大 $\tan\beta$ Higgs 介导贡献尤其敏感，而 muon $g-2$ 又会对较轻的电弱超对称粒子施加互补要求。

对于 NMSSM，情况与 MSSM 类似，但由于额外 singlet 标量和赝标量的存在，Higgs sector 的 flavour 影响可能更加丰富，因此 $B_s\to\mu^+\mu^-$、$B\to X_s\gamma$ 和若干半轻子 observables 都可能给出额外限制。对于模型无关 EFT 分析，[modelindep_chi2.c](modelindep_chi2.c) 直接把 $\Delta C_9^e$、$\Delta C_9^\mu$、$\Delta C_{10}^e$ 和 $\Delta C_{10}^\mu$ 当作输入，这样用户无需先指定完整的 UV 模型，就可以直接研究实验数据更偏好哪一种有效算符偏移。这一做法在 flavour anomaly 文献中非常常见，因为它能够把“数据偏好什么算符结构”和“哪个具体模型能产生这种结构”两个问题分开处理。

## 十一、统计意义上的“限制”在程序里是怎样量化的

在日常讨论中，人们常说某个 flavour observable “限制了” 某个模型参数，但在代码中，这种限制并不是一句定性话，而是通过具体统计量定义出来的。SuperIso 通过 [src/chi2.c](src/chi2.c) 构造理论与实验的总协方差矩阵，再计算全局 $\chi^2$。如果某个参数点给出的预测与实验中心值普遍接近，且落在协方差所允许的误差范围之内，那么该点的 $\chi^2$ 就较小，说明它与数据相容。相反，如果某个参数点同时在多个观测量上偏离实验，或者在关键观测量上偏离得非常明显，那么全局 $\chi^2$ 就会迅速变大，从而说明该点受到数据排斥。

这一点尤其重要，因为 flavour observable 之间不是彼此独立的。很多观测量共享同样的 CKM 参数、相同的形状因子，或者来自同一个实验分析，因此它们的误差存在相关性。若忽略这些相关性，参数限制就可能过紧或过松。SuperIso 通过 nuisance 参数和协方差矩阵来处理这一问题，使得得到的限制更接近严格的全局拟合结果，而不是简单的单观测量卡方叠加。

## 十二、项目在实际研究中的意义

综合来看，SuperIso 的价值在于它提供了一条完整而清晰的研究路径。研究者可以从一个具体模型出发，通过外部谱生成器或用户提供的 SLHA 文件获得低能参数；然后利用 SuperIso 计算 flavour observable 和精密量；再利用这些观测量与实验比较，判断该模型参数点是否仍然可行。如果需要更进一步的系统分析，还可以通过 [sm_chi2.c](sm_chi2.c)、[slha_chi2.c](slha_chi2.c) 或 [modelindep_chi2.c](modelindep_chi2.c) 做标准模型、具体新物理点或模型无关 Wilson 系数的全局 $\chi^2$ 拟合。

从代码角度看，这个项目的结构是相当成熟的：顶层主程序负责场景入口，`src/` 中的物理模块负责具体观测量计算，`interpreter` 负责统一 observable 接口，`chi2_input` 负责实验与 nuisance 数据，`chi2` 模块负责统计整合。从物理角度看，它覆盖了 flavour phenomenology 中最重要的一批观测量，尤其适合研究 $b\to s\gamma$、$b\to s\ell\ell$、$B_s\to\mu^+\mu^-$、$B\to D^{(*)}\tau\nu$、稀有 kaon 衰变以及 muon $g-2$ 对 2HDM、MSSM、NMSSM 和模型无关 EFT 的限制。因此，如果要把这个项目用一句话概括，可以说：它是一个把高能模型参数系统地投影到低能 flavour observable，并通过全局统计方法给出新物理限制的专业工具。