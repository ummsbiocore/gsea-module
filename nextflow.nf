$HOSTNAME = ""
params.outdir = 'results'  


if (!params.run_gsea){params.run_gsea = ""} 
if (!params.input){params.input = ""} 
if (!params.postfix){params.postfix = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)

Channel.value(params.run_gsea).set{g_1_2_g_0}
g_2_1_g_0 = file(params.input, type: 'any')
Channel.value(params.postfix).set{g_5_0_g_0}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 30 
}
//* platform
//* platform
//* autofill

process Prepare_GSEA {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /GSEA$/) "GSEA/$filename"}
input:
 val postfix
 path input
 val run_GSEA

output:
 path "GSEA_reports"  ,emit:g_0_outputFile00 
 path "GSEA"  ,emit:g_0_outputFile11 

container 'quay.io/viascientific/gsea_module:1.0.0'

// SET second output to "{GSEA,outputs}" when launched apps can reach parent directory

when:
run_GSEA == 'yes'

script:

if (postfix.equals(" ")) {
	postfix = ''
}

event_column = params.Prepare_GSEA.event_column
fold_change_column = params.Prepare_GSEA.fold_change_column

if (params.genome_build.startsWith("human_")) {
	species = 'human'
} else if (params.genome_build.startsWith("mouse_")) {
	species = 'mouse'
} else {
	species = 'NA'
}

local_species = params.Prepare_GSEA.local_species
H = params.Prepare_GSEA.H
C1 = params.Prepare_GSEA.C1
C2 = params.Prepare_GSEA.C2
C2_CGP = params.Prepare_GSEA.C2_CGP
C2_CP = params.Prepare_GSEA.C2_CP
C3_MIR = params.Prepare_GSEA.C3_MIR
C3_TFT = params.Prepare_GSEA.C3_TFT
C4 = params.Prepare_GSEA.C4
C4_CGN = params.Prepare_GSEA.C4_CGN
C4_CM = params.Prepare_GSEA.C4_CM
C5 = params.Prepare_GSEA.C5
C5_GO = params.Prepare_GSEA.C5_GO
C5_GO_BP = params.Prepare_GSEA.C5_GO_BP
C5_GO_CC = params.Prepare_GSEA.C5_GO_CC
C5_GO_MF = params.Prepare_GSEA.C5_GO_MF
C5_HPO = params.Prepare_GSEA.C5_HPO
C6 = params.Prepare_GSEA.C6
C7 = params.Prepare_GSEA.C7
C8 = params.Prepare_GSEA.C8
MH = params.Prepare_GSEA.MH
M1 = params.Prepare_GSEA.M1
M2 = params.Prepare_GSEA.M2
M2_CGP = params.Prepare_GSEA.M2_CGP
M2_CP = params.Prepare_GSEA.M2_CP
M3_GTRD = params.Prepare_GSEA.M3_GTRD
M3_miRDB = params.Prepare_GSEA.M3_miRDB
M5 = params.Prepare_GSEA.M5
M5_GO = params.Prepare_GSEA.M5_GO
M5_GO_BP = params.Prepare_GSEA.M5_GO_BP
M5_GO_CC = params.Prepare_GSEA.M5_GO_CC
M5_GO_MF = params.Prepare_GSEA.M5_GO_MF
M5_MPT = params.Prepare_GSEA.M5_MPT
M8 = params.Prepare_GSEA.M8

minSize = params.Prepare_GSEA.minSize
maxSize = params.Prepare_GSEA.maxSize

nes_significance_cutoff = params.Prepare_GSEA.nes_significance_cutoff
padj_significance_cutoff = params.Prepare_GSEA.padj_significance_cutoff

seed = params.Prepare_GSEA.seed

//* @style @condition:{local_species="human",H,C1,C2,C2_CGP,C2_CP,C3_MIR,C3_TFT,C4,C4_CGN,C4_CM,C5,C5_GO,C5_GO_BP,C5_GO_CC,C5_GO_MF,C5_HPO,C6,C7,C8},{local_species="mouse",MH,M1,M2,M2_CGP,M2_CP,M3_GTRD,M3_miRDB,M5,M5_GO,M5_GO_BP,M5_GO_CC,M5_GO_MF,M5_MPT,M8} @multicolumn:{event_column, fold_change_column}, {gmt_list, minSize, maxSize},{H,C1,C2,C2_CGP}, {C2_CP,C3_MIR,C3_TFT,C4},{C4_CGN,C4_CM, C5,C5_GO},{C5_GO_BP,C5_GO_CC,C5_GO_MF,C5_HPO},{C6,C7,C8},{MH,M1,M2,M2_CGP},{M2_CP,M3_GTRD,M3_miRDB,M5},{M5_GO,M5_GO_BP,M5_GO_CC,M5_GO_MF},{M5_MPT,M8},{nes_significance_cutoff, padj_significance_cutoff}

H        = H        == 'true' && species == 'human' ? ' h.all.v2023.1.Hs.entrez.gmt'    : ''
C1       = C1       == 'true' && species == 'human' ? ' c1.all.v2023.1.Hs.entrez.gmt'   : ''
C2       = C2       == 'true' && species == 'human' ? ' c2.all.v2023.1.Hs.entrez.gmt'   : ''
C2_CGP   = C2_CGP   == 'true' && species == 'human' ? ' c2.cgp.v2023.1.Hs.entrez.gmt'   : ''
C2_CP    = C2_CP    == 'true' && species == 'human' ? ' c2.cp.v2023.1.Hs.entrez.gmt'    : ''
C3_MIR   = C3_MIR   == 'true' && species == 'human' ? ' c3.mir.v2023.1.Hs.entrez.gmt'   : ''
C3_TFT   = C3_TFT   == 'true' && species == 'human' ? ' c3.tft.v2023.1.Hs.entrez.gmt'   : ''
C4       = C4       == 'true' && species == 'human' ? ' c4.all.v2023.1.Hs.entrez.gmt'   : ''
C4_CGN   = C4_CGN   == 'true' && species == 'human' ? ' c4.cgn.v2023.1.Hs.entrez.gmt'   : ''
C4_CM    = C4_CM    == 'true' && species == 'human' ? ' c4.cm.v2023.1.Hs.entrez.gmt'    : ''
C5       = C5       == 'true' && species == 'human' ? ' c5.all.v2023.1.Hs.entrez.gmt'   : ''
C5_GO    = C5_GO    == 'true' && species == 'human' ? ' c5.go.v2023.1.Hs.entrez.gmt'    : ''
C5_GO_BP = C5_GO_BP == 'true' && species == 'human' ? ' c5.go.bp.v2023.1.Hs.entrez.gmt' : ''
C5_GO_CC = C5_GO_CC == 'true' && species == 'human' ? ' c5.go.cc.v2023.1.Hs.entrez.gmt' : ''
C5_GO_MF = C5_GO_MF == 'true' && species == 'human' ? ' c5.go.mf.v2023.1.Hs.entrez.gmt' : ''
C5_HPO   = C5_HPO   == 'true' && species == 'human' ? ' c5.hpo.v2023.1.Hs.entrez.gmt'   : ''
C6       = C6       == 'true' && species == 'human' ? ' c6.all.v2023.1.Hs.entrez.gmt'   : ''
C7       = C7       == 'true' && species == 'human' ? ' c7.all.v2023.1.Hs.entrez.gmt'   : ''
C8       = C8       == 'true' && species == 'human' ? ' c8.all.v2023.1.Hs.entrez.gmt'   : ''
MH       = MH       == 'true' && species == 'mouse' ? ' mh.all.v2023.1.Mm.entrez.gmt'   : ''    
M1       = M1       == 'true' && species == 'mouse' ? ' m1.all.v2023.1.Mm.entrez.gmt'   : ''    
M2       = M2       == 'true' && species == 'mouse' ? ' m2.all.v2023.1.Mm.entrez.gmt'   : ''    
M2_CGP   = M2_CGP   == 'true' && species == 'mouse' ? ' m2.cgp.v2023.1.Mm.entrez.gmt'   : ''
M2_CP    = M2_CP    == 'true' && species == 'mouse' ? ' m2.cp.v2023.1.Mm.entrez.gmt'    : '' 
M3_GTRD  = M3_GTRD  == 'true' && species == 'mouse' ? ' m3.gtrd.v2023.1.Mm.entrez.gmt'  : ''
M3_miRDB = M3_miRDB == 'true' && species == 'mouse' ? ' m3.mirdb.v2023.1.Mm.entrez.gmt' : ''
M5       = M5       == 'true' && species == 'mouse' ? ' m5.all.v2023.1.Mm.entrez.gmt'   : ''    
M5_GO    = M5_GO    == 'true' && species == 'mouse' ? ' m5.go.v2023.1.Mm.entrez.gmt'    : '' 
M5_GO_BP = M5_GO_BP == 'true' && species == 'mouse' ? ' m5.go.bp.v2023.1.Mm.entrez.gmt' : ''
M5_GO_CC = M5_GO_CC == 'true' && species == 'mouse' ? ' m5.go.cc.v2023.1.Mm.entrez.gmt' : ''
M5_GO_MF = M5_GO_MF == 'true' && species == 'mouse' ? ' m5.go.mf.v2023.1.Mm.entrez.gmt' : ''
M5_MPT   = M5_MPT   == 'true' && species == 'mouse' ? ' m5.mpt.v2023.1.Mm.entrez.gmt'   : ''
M8       = M8       == 'true' && species == 'mouse' ? ' m8.all.v2023.1.Mm.entrez.gmt'   : ''

gmt_list = H + C1 + C2 + C2_CGP + C2_CP + C3_MIR + C3_TFT + C4 + C4_CGN + C4_CM + C5 + C5_GO + C5_GO_BP + C5_GO_CC + C5_GO_MF + C5_HPO + C6 + C7 + C8 + MH + M1 + M2 + M2_CGP + M2_CP + M3_GTRD + M3_miRDB + M5 + M5_GO + M5_GO_BP + M5_GO_CC + M5_GO_MF + M5_MPT + M8

if (gmt_list.equals("")){
	gmt_list = 'h.all.v2023.1.Hs.entrez.gmt'
}

"""
prepare_GSEA.py \
--input ${input} --species ${species} --event-column ${event_column} --fold-change-column ${fold_change_column} \
--GMT-key /data/gmt_key.txt --GMT-source /data --GMT-list ${gmt_list} --minSize ${minSize} --maxSize ${maxSize} \
--NES ${nes_significance_cutoff} --pvalue ${padj_significance_cutoff} \
--seed ${seed} --threads ${task.cpus} --postfix '_gsea${postfix}'

cp -R outputs GSEA/

mkdir GSEA_reports
cp -R GSEA GSEA_reports/
"""

}


workflow {

g_5_0_g_0= g_5_0_g_0.ifEmpty("") 


Prepare_GSEA(g_5_0_g_0,g_2_1_g_0,g_1_2_g_0)
g_0_outputFile00 = Prepare_GSEA.out.g_0_outputFile00
g_0_outputFile11 = Prepare_GSEA.out.g_0_outputFile11


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
