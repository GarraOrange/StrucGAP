from pathlib import Path
import os
import sys
# ============================================================
# pyecharts 本地资源设置：解决 PhantomJS 找不到 echarts 的问题
# ============================================================
import pyecharts
from pyecharts.globals import CurrentConfig

PYECHARTS_ASSET_DIR = Path(r"D:\Research\tools\pyecharts-assets")
CurrentConfig.ONLINE_HOST = PYECHARTS_ASSET_DIR.resolve().as_uri() + "/"

print("[pyecharts version]", pyecharts.__version__)
print("[pyecharts asset dir]", CurrentConfig.ONLINE_HOST)
from strucgap.preprocess import StrucGAP_Preprocess
from strucgap.glycanstructure import StrucGAP_GlycanStructure
from strucgap.glycosite import StrucGAP_GlycoSite
from strucgap.glycopeptidequant import StrucGAP_GlycoPeptideQuant
from strucgap.functionannotation import StrucGAP_FunctionAnnotation
from strucgap.datavisualization import StrucGAP_DataVisualization
from strucgap.insighttracker import StrucGAP_InsightTracker

# 0. 路径和参数：一般只需要改这里
BASE_DIR = Path(r"D:\Research\lvy\StrucGAP\test")
OUT_DIR = BASE_DIR / "output"

PGLYCO_STRUCTURE_FILE = BASE_DIR / "pd structure.xlsx"
QUANT_FILE = BASE_DIR / "s6.xlsx"
SAMPLE_GROUP_FILE = BASE_DIR / "sample_group.xlsx"


BRANCH_FILE_CANDIDATES = list(BASE_DIR.glob("branch_structures*.xlsx"))
BRANCH_FILE = BRANCH_FILE_CANDIDATES[0] if BRANCH_FILE_CANDIDATES else BASE_DIR / "branch_structures_18_mice uterus.0240401.xlsx"
WURCS_FILE = BASE_DIR / "glycosmos_glycans_wurcs.csv"

DATA_SHEET_NAME = "Sheet1"
QUANT_SHEET_NAME = "Sheet1"

# 官方 pGlyco3 示例定量列名
QUANT_COLS = [
    "Young-SN-1", "Young-SN-2", "Young-SN-3",
    "PD-SN-1", "PD-SN-2", "PD-SN-3",
]

# 官方示例中的 abundance_ratio。
# 注意：这里有 10 个数，是官方示例参数
ABUNDANCE_RATIO = [
    1.240003449, 0, 1.344387558, 0, 1.576533442,
    0, 1, 0, 1.956346409, 1.517000766,
]

# 是否运行两个流程
RUN_STRUCTURE_ONLY = False
RUN_STRUCTURE_WITH_QUANT = True

# 差异分析阈值：官方示例使用 fc=4.2
DIFF_FC = 4.2



# ============================================================
# 1. 工具函数
# ============================================================

def check_files() -> None:
    """检查必要输入文件是否存在。"""
    required_files = {
        "pGlyco3 structure file": PGLYCO_STRUCTURE_FILE,
        "quantification file": QUANT_FILE,
        "sample group file": SAMPLE_GROUP_FILE,
        "branch list file": BRANCH_FILE,
        "GlyTouCan WURCS file": WURCS_FILE,
    }

    missing = []
    for name, path in required_files.items():
        if not path.exists():
            missing.append(f"{name}: {path}")

    if missing:
        print("\n[ERROR] 缺少必要文件：")
        for item in missing:
            print("  -", item)
        print("\n请把缺失文件放到：", BASE_DIR)
        print("其中 branch list 通常来自官方 tests 目录：branch_structures_18_mice uterus.0240401.xlsx")
        print("GlyTouCan WURCS 文件通常来自官方 tests 目录：glycosmos_glycans_wurcs.csv")
        sys.exit(1)


def make_and_enter_output_dir(subdir: str) -> Path:
    """创建并进入某个输出子目录。StrucGAP 的 output() 通常写入当前工作目录。"""
    work_dir = OUT_DIR / subdir
    work_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(work_dir)
    print(f"\n[Working directory] {work_dir}")
    return work_dir


def run_structure_only(data_manager: StrucGAP_InsightTracker) -> None:
    """
    只使用 pGlyco3 结构信息：
    Preprocess -> GlycanStructure -> GlycoSite
    """
    make_and_enter_output_dir("01_structure_only")

    print("\n========== 01_structure_only: Preprocess ==========")
    module1 = StrucGAP_Preprocess(
        data_dir=str(PGLYCO_STRUCTURE_FILE),
        data_sheet_name=DATA_SHEET_NAME,
        sample_group_data_dir=str(SAMPLE_GROUP_FILE),
        branch_list_dir=str(BRANCH_FILE),
        data_manager=data_manager,
        search_engine="pGlyco3",
    )

    module1.data_cleaning(data_type="tmt")
    module1.cv_raw(threshold="no", fc_recommendation=False)
    module1.fdr(feature_type="no")
    module1.outliers(abundance_ratio=ABUNDANCE_RATIO)
    module1.cv(threshold="no")
    module1.psm(psm_number="no", fc_recommendation=False)
    module1.annotation(
        glytoucan=True,
        glytoucan_structure=True,
        glytoucan_wurcs_file=str(WURCS_FILE),
        biosynthetic_pathways=True,
        glycobiology_filter=True,
    )
    module1.output()

    print("\n========== 01_structure_only: GlycanStructure ==========")
    module2 = StrucGAP_GlycanStructure(
        gs_data=module1,
        data_manager=data_manager,
        data_type="psm_filtered",
    )
    module2.statistics(remove_oligo_mannose=False)
    module2.structure_statistics()
    module2.lacdinac()
    module2.cor()
    module2.isoforms()
    module2.output()

    print("\n========== 01_structure_only: GlycoSite ==========")
    module3 = StrucGAP_GlycoSite(module1, data_manager=data_manager)
    module3.glycoprotein_site()
    module3.glycopeptide_site()
    module3.specific_site()
    module3.output()

    print("\n[OK] structure-only workflow finished.")


def run_structure_with_quant(data_manager: StrucGAP_InsightTracker):
    """
    使用 pGlyco3 结构信息 + 外部定量表：
    Preprocess -> GlycanStructure -> GlycoSite -> GlycoPeptideQuant -> FunctionAnnotation
    """
    make_and_enter_output_dir("02_structure_with_quant")

    print("\n========== 02_structure_with_quant: Preprocess ==========")
    module1 = StrucGAP_Preprocess(
        data_dir=str(PGLYCO_STRUCTURE_FILE),
        data_sheet_name=DATA_SHEET_NAME,
        sample_group_data_dir=str(SAMPLE_GROUP_FILE),
        branch_list_dir=str(BRANCH_FILE),
        data_manager=data_manager,
        search_engine="pGlyco3",
    )

    module1.data_cleaning(
        data_type="tmt",
        quantification_from_no_strucgp=True,
        quantification_data_dir=str(QUANT_FILE),
        sheet_name=QUANT_SHEET_NAME,
        quant_cols=QUANT_COLS,
    )
    module1.cv_raw(threshold="no", fc_recommendation=True)
    module1.fdr(feature_type="no")
    module1.outliers(abundance_ratio=ABUNDANCE_RATIO)
    module1.cv(threshold="no")
    module1.psm(psm_number="no", fc_recommendation=True)
    module1.annotation(
        glytoucan=True,
        glytoucan_structure=True,
        glytoucan_wurcs_file=str(WURCS_FILE),
        biosynthetic_pathways=True,
        glycobiology_filter=True,
    )
    module1.output()

    print("\n========== 02_structure_with_quant: GlycanStructure ==========")
    module2 = StrucGAP_GlycanStructure(
        gs_data=module1,
        data_manager=data_manager,
        data_type="psm_filtered",
    )
    module2.statistics(remove_oligo_mannose=False)
    module2.structure_statistics()
    module2.lacdinac()
    module2.cor()
    module2.isoforms()
    module2.output()

    print("\n========== 02_structure_with_quant: GlycoSite ==========")
    module3 = StrucGAP_GlycoSite(module1, data_manager=data_manager)
    module3.glycoprotein_site()
    module3.glycopeptide_site()
    module3.specific_site()
    module3.output()

    print("\n========== 02_structure_with_quant: GlycoPeptideQuant ==========")
    module4 = StrucGAP_GlycoPeptideQuant(
        module1,
        data_type="psm_filtered",
        data_manager=data_manager,
    )
    module4.statistics()
    module4.statistics_index()
    module4.differential_analysis(pvalue_type="pvalue_ttest", fc=DIFF_FC)
    module4.threshold_variation_analysis(
        pvalue_type="pvalue_ttest",
        statistic_index="fc",
        fc_range=[3, 6, 10, 50, 100],
    )
    module4.glycopeptide_glycosite_glycan_variation(fc=DIFF_FC)
    module4.glycoprotein_glycosite_glycan_variation(fc=DIFF_FC)
    module4.output()

    data_manager.key_information_extraction(module="StrucGAP_GlycoPeptideQuant")

    print("\n========== 02_structure_with_quant: FunctionAnnotation based on module1 ==========")
    module5 = StrucGAP_FunctionAnnotation(module1, data_manager=data_manager)
    module5.ora(
        organism="mmusculus",
        background_input=False,
        up_down_fc_threshold=1.5,
    )
    module5.go_function_structure(function_data="ora_no_background_both_proteins_result")
    module5.output()
    data_manager.key_information_extraction(module="StrucGAP_FunctionAnnotation")

    print("\n========== 02_structure_with_quant: FunctionAnnotation based on module4 ==========")
    module5 = StrucGAP_FunctionAnnotation(module4, data_manager=data_manager)
    module5.ora(
        organism="mmusculus",
        selected_terms=["GO:MF", "GO:CC", "GO:BP"],
        enrich_feature="glycopeptide",
        background_input=False,
        up_down_fc_threshold=DIFF_FC,
        pvalue_type="pvalue_ttest",
    )

    module5.go_function_structure(function_data="ora_no_background_up_result")
    module5.output()
    data_manager.key_information_extraction(module="StrucGAP_FunctionAnnotation")

    module5.go_function_structure(function_data="ora_no_background_down_result")
    module5.output()
    data_manager.key_information_extraction(module="StrucGAP_FunctionAnnotation")

    print("\n[OK] structure-with-quant workflow finished.")
    return module1, module2, module3, module4


# ============================================================
# 2. 主流程
# ============================================================

def main() -> None:
    print("[Base directory]", BASE_DIR)
    print("[Output directory]", OUT_DIR)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    check_files()

    data_manager = StrucGAP_InsightTracker()

    module1 = None
    module2 = None
    module3 = None
    module4 = None

    if RUN_STRUCTURE_ONLY:
        run_structure_only(data_manager)

    if RUN_STRUCTURE_WITH_QUANT:
        module1, module2, module3, module4 = run_structure_with_quant(data_manager)

    # 保存 data_manager 状态
    make_and_enter_output_dir("99_data_manager")
    data_manager.output_pickle()

    # 官方可视化脚本需要 module1/module3/module4/module7
    if module1 is None or module3 is None or module4 is None:
        raise RuntimeError(
            "官方可视化需要先运行 RUN_STRUCTURE_WITH_QUANT=True，"
            "因为 analysis_report_pglyco3_local.py 依赖 module1/module3/module4。"
        )

    module7 = StrucGAP_DataVisualization(data_manager=data_manager)

    report_script = BASE_DIR / "analysis_report_pglyco3_local.py"

    exec_context = globals().copy()
    exec_context.update({
        "module1": module1,
        "module2": module2,
        "module3": module3,
        "module4": module4,
        "module7": module7,
        "data_manager": data_manager,
    })

    exec(report_script.read_text(encoding="utf-8"), exec_context)

    print("\nAll workflows finished.")
    print("Results are under:", OUT_DIR)


if __name__ == "__main__":
    main()
