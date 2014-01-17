#!/bin/sh

samples="pY2_EV_exp_nucleosome pY2_EV_stat_nucleosome pY2_Mm3a_exp_nucleosome pY2_Mm3a_stat_nucleosome pY2_Mm3b_exp_nucleosome pY2_Mm3b_stat_nucleosome pY2_EV_exp_naked pY2_Mm3a_exp_naked pY2_Mm3b_exp_naked pY2_Mm3b_stat_naked"

### sam2bed
sam2bed() {
cd /u/home/k/kaifu/project-mcdb/nuc_methy/rawdata_nuc/sam
for each in $samples;
do
python ../../phase1_nuc_analysis/bin/sam2bed.py $each.sam ../bed/$each.bed
done
}

### get unique reads
uniquebed() {
cd /u/home/k/kaifu/project-mcdb/nuc_methy/rawdata_nuc/bed
for each in $samples;
do
python /u/home/k/kaifu/project-mcdb/nuc_methy/phase1_nuc_analysis/bin/get_uniq_read.py -i $each.bed -o ../uniquebed/$each.uniq.bed
done
}

### average profile
ave_profile() {
cd /u/home/k/kaifu/project-mcdb/nuc_methy/phase1_nuc_analysis/ave_profile
gene=/u/home/k/kaifu/project-mcdb/nuc_methy/genome_infor/sgd_gene.bed

python ../bin/nucleosome_profile_all.py -i ../../rawdata_nuc/uniquebed/pY2_EV_exp_nucleosome_uniq.bed --tss EV_exp_nucleosome_tss --tts EV_exp_nucleosome_tts --tsseach EV_exp_nucleosome_tss_each --ttseach EV_exp_nucleosome_tts_each -s yeast -g $gene --name EV_exp_nucleosome
python ../bin/nucleosome_profile_all.py -i ../../rawdata_nuc/uniquebed/pY2_EV_stat_nucleosome_uniq.bed --tss EV_stat_nucleosome_tss --tts EV_stat_nucleosome_tts --tsseach EV_stat_nucleosome_tss_each --ttseach EV_stat_nucleosome_tts_each -s yeast -g $gene --name EV_stat_nucleosome
python ../bin/nucleosome_profile_all.py -i ../../rawdata_nuc/uniquebed/pY2_Mm3a_exp_nucleosome_uniq.bed --tss Mm3a_exp_nucleosome_tss --tts Mm3a_exp_nucleosome_tts --tsseach Mm3a_exp_nucleosome_tss_each --ttseach Mm3a_exp_nucleosome_tts_each -s yeast -g $gene --name Mm3a_exp_nucleosome
python ../bin/nucleosome_profile_all.py -i ../../rawdata_nuc/uniquebed/pY2_Mm3a_stat_nucleosome_uniq.bed --tss Mm3a_stat_nucleosome_tss --tts Mm3a_stat_nucleosome_tts --tsseach Mm3a_exp_nucleosome_tss_each --ttseach Mm3a_exp_nucleosome_tts_each -s yeast -g $gene --name Mm3a_exp_nucleosome
python ../bin/nucleosome_profile_all.py -i ../../rawdata_nuc/uniquebed/pY2_Mm3b_exp_nucleosome_uniq.bed --tss Mm3b_exp_nucleosome_tss --tts Mm3b_exp_nucleosome_tts --tsseach Mm3b_exp_nucleosome_tss_each --ttseach Mm3b_exp_nucleosome_tts_each -s yeast -g $gene --name Mm3b_exp_nucleosome
python ../bin/nucleosome_profile_all.py -i ../../rawdata_nuc/uniquebed/pY2_Mm3b_stat_nucleosome_uniq.bed --tss Mm3b_stat_nucleosome_tss --tts Mm3b_stat_nucleosome_tts --tsseach Mm3b_stat_nucleosome_tss_each --ttseach Mm3b_stat_nucleosome_tts_each -s yeast -g $gene --name Mm3b_stat_nucleosome
}

### 10bp periodicity
dinucleotide() {
cd /u/home/k/kaifu/project-mcdb/nuc_methy/phase1_nuc_analysis/10bp_periodicity
for each in $samples;
do 
python /u/home/k/kaifu/project-mcdb/nuc_methy/phase1_nuc_analysis/bin/Nuc_code/Nuc_code_detection.py ../../rawdata_nuc/uniquebed/$each.uniq.bed -s yeast -o $each.nuc_code
done
}

### statistical positioning
stat_position() {
cd /u/home/k/kaifu/project-mcdb/nuc_methy/phase1_nuc_analysis/statistical_positioning
for each in $samples;
do
python ../bin/Statistical_Positioning/Nuc_positioning_relationship.py ../../rawdata_nuc/uniquebed/$each.uniq.bed ../../rawdata_nuc/uniquebed/$each.uniq.bed -s yeast -o $each.stat_position
done
}

#sam2bed
#uniquebed
#ave_profile
dinucleotide
stat_position


