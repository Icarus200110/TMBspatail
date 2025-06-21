#import pysam

def get_tgp_frequence(chromo, position, frequence_type):
    gvcf = "tgp_phase3_small_vars.vcf.gz"
    vcf_in = pysam.VariantFile(gvcf)
    chromo = chromo[3:]
    frequence_type = frequence_type + '='
    try:
        for rec in vcf_in.fetch(chromo, position - 1, position):
            this_base = str(rec)
            break
        this_base_array = this_base.split()
        ref_length = len(this_base_array[3])
        alt = this_base_array[4]
        this_base_info_array = this_base_array[7].split(';')
        for this_base_info_element in this_base_info_array:
            if this_base_info_element.startswith(frequence_type):
                frequence = this_base_info_element.split('=')[-1]
                return ref_length, alt, frequence
    except:
        return ""
