from vannotplus.family.barcode import main_barcode
from vannotplus.exomiser.exomiser import main_exomiser
from vannotplus.annot.score import main_annot
from vannotplus.__main__ import load_config, main_config

vcf_file = "/STARK/output/tmp/tmp_230908_NB551027_0178_AHW7KLBGXT/SGT2301624.vcf.gz"
output = "/STARK/output/tmp/tmp_230908_NB551027_0178_AHW7KLBGXT/output.vcf.gz"
config = "/STARK/config/variantannotation/vannot/vannotplus.yml"
main_exomiser(vcf_file, output, "WES_AGILENT", load_config(config))
