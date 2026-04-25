from pathlib import Path
import os
import time

from strucgap.preprocess import StrucGAP_Preprocess
from strucgap.glycopeptidequant import StrucGAP_GlycoPeptideQuant
from strucgap.functionannotation import StrucGAP_FunctionAnnotation
from strucgap.insighttracker import StrucGAP_InsightTracker


#目录
BASE_DIR = Path(r"D:\Research\lvy\StrucGAP\test")
OUT_DIR = BASE_DIR / "output"
WORK_DIR = OUT_DIR / "02_structure_with_quant"

PGLYCO_STRUCTURE_FILE = BASE_DIR / "pd structure.xlsx"
QUANT_FILE = BASE_DIR / "s6.xlsx"
SAMPLE_GROUP_FILE = BASE_DIR / "sample_group.xlsx"
WURCS_FILE = BASE_DIR / "glycosmos_glycans_wurcs.csv"

BRANCH_FILE_CANDIDATES = list(BASE_DIR.glob("branch_structures*.xlsx"))
if len(BRANCH_FILE_CANDIDATES) == 0:
    raise FileNotFoundError("没有找到 branch_structures*.xlsx，请确认 branch 文件在 BASE_DIR 下。")
BRANCH_FILE = BRANCH_FILE_CANDIDATES[0]


# ==================================================
# 1. 参数
# ==================================================
DATA_SHEET_NAME = "Sheet1"
QUANT_SHEET_NAME = "Sheet1"

QUANT_COLS = [
    "Young-SN-1",
    "Young-SN-2",
    "Young-SN-3",
    "PD-SN-1",
    "PD-SN-2",
    "PD-SN-3",
]

ABUNDANCE_RATIO = [
    1.240003449, 0, 1.344387558, 0, 1.576533442,
    0, 1, 0, 1.956346409, 1.517000766,
]

FC_CUTOFF = 4.2
ORGANISM = "mmusculus"


# ==================================================
# 2. 小工具：带重试的 g:Profiler
# ==================================================
def run_with_retry(func, max_retry=3, sleep_seconds=20):
    last_error = None

    for i in range(max_retry):
        try:
            print(f"[Try {i + 1}/{max_retry}] running online ORA ...")
            return func()
        except Exception as e:
            last_error = e
            print(f"[Warning] failed at attempt {i + 1}: {e}")
            if i < max_retry - 1:
                print(f"[Wait] sleep {sleep_seconds} seconds and retry ...")
                time.sleep(sleep_seconds)

    raise last_error


# ==================================================
# 3. 最小重建 module1
# ==================================================
def rebuild_module1(data_manager):
    print("\n========== Rebuild minimal module1: Preprocess ==========")

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

    return module1


# ==================================================
# 4. 最小重建 module4
# ==================================================
def rebuild_module4(module1, data_manager):
    print("\n========== Rebuild minimal module4: GlycoPeptideQuant ==========")

    module4 = StrucGAP_GlycoPeptideQuant(
        module1,
        data_type="psm_filtered",
        data_manager=data_manager,
    )

    module4.statistics()
    module4.statistics_index()
    module4.differential_analysis(
        pvalue_type="pvalue_ttest",
        fc=FC_CUTOFF,
    )

    # 注意：下面这些都先不跑，因为它们已经在上一次输出过，
    # 且对继续 FunctionAnnotation 不是必须。
    #
    # module4.threshold_variation_analysis(...)
    # module4.glycopeptide_glycosite_glycan_variation(...)
    # module4.glycoprotein_glycosite_glycan_variation(...)
    # module4.output()

    return module4


# ==================================================
# 5.FunctionAnnotation：glycopeptide up/down
# ==================================================
def continue_function_annotation(module4, data_manager):
    print("\n========== Continue module5: FunctionAnnotation ==========")

    module5 = StrucGAP_FunctionAnnotation(
        module4,
        data_manager=data_manager,
    )

    # 5.1 ORA
    def do_ora():
        return module5.ora(
            organism=ORGANISM,
            selected_terms=["GO:MF", "GO:CC", "GO:BP"],
            enrich_feature="glycopeptide",
            background_input=False,
            up_down_fc_threshold=FC_CUTOFF,
            pvalue_type="pvalue_ttest",
        )

    run_with_retry(do_ora, max_retry=3, sleep_seconds=20)

    # 5.2 up result
    print("\n========== GO structure: up result ==========")
    module5.go_function_structure(function_data="ora_no_background_up_result")
    module5.output()

    data_manager.key_information_extraction(
        module="StrucGAP_FunctionAnnotation"
    )

    # 5.3 down result
    print("\n========== GO structure: down result ==========")
    module5.go_function_structure(function_data="ora_no_background_down_result")
    module5.output()

    data_manager.key_information_extraction(
        module="StrucGAP_FunctionAnnotation"
    )

    return module5


# ==================================================
# 6. 主函数
# ==================================================
def main():
    print(f"[BASE_DIR] {BASE_DIR}")
    print(f"[WORK_DIR] {WORK_DIR}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    os.chdir(WORK_DIR)

    data_manager = StrucGAP_InsightTracker()

    module1 = rebuild_module1(data_manager)
    module4 = rebuild_module4(module1, data_manager)
    module5 = continue_function_annotation(module4, data_manager)

    print("\n========== Save resume data_manager ==========")
    data_manager.output_pickle()

    print("\n[DONE] Resume FunctionAnnotation finished.")


if __name__ == "__main__":
    main()