#!/bin/bash
#SBATCH --job-name=vers2.bw
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long




mapping=/home/ualtintas/LncapExo/CAGE


outData=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/Cage/bigWig



srx=(
"SRX882919"
"SRX882920"
"SRX882921"
"SRX882922"
"SRX882923"
"SRX882924"
"SRX882925"
"SRX882926"
"SRX882927"
"SRX882928"
"SRX882929"
"SRX882930"
"SRX882931"
"SRX882932"
)


merge=(
"LNCaP.0h.rep1"
"LNCaP.3h.rep1"
"LNCaP.6h.rep1"
"LNCaP.12h.rep1"
"LNCaP.18h.rep1"
"LNCaP.48h.rep1"
"LNCaP.72h.rep1"
)



for i in {0..6}
do
  SRX1=${srx[$i]};
  SRX2=${srx[$i+7]};
  merg=${merge[$i]};


  bamSrx1=$mapping/$SRX1".processed.bam";
  bamSrx2=$mapping/$SRX2".processed.bam";
  outMerg=$mapping/$merg".bam";
  outSrxP=$outData/$merg".-.BPM.bw";
  outSrxM=$outData/$merg".+.BPM.bw";

  echo $merg $i $bamSrx1 $bamSrx2 ;

  samtools merge -f -@ 10 $outMerg $bamSrx1 $bamSrx2;
  samtools index $outMerg;

  bamCoverage --bam $outMerg -o $outSrxP --normalizeUsing BPM -p 30 --filterRNAstrand forward -bs 20
  bamCoverage --bam $outMerg -o $outSrxM --normalizeUsing BPM -p 30 --filterRNAstrand reverse -bs 20
done



samtools merge -f -@ 10 $mapping/cageMerged.bam $mapping/*.processed.bam
samtools index $mapping/cageMerged.bam

bamCoverage --bam $mapping/cageMerged.bam -o $outData/cageMerged.-.BPM.bw --normalizeUsing BPM -p 30 --filterRNAstrand forward -bs 20
bamCoverage --bam $mapping/cageMerged.bam -o $outData/cageMerged.+.BPM.bw --normalizeUsing BPM -p 30 --filterRNAstrand reverse -bs 20
