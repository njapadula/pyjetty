#!/bin/bash

cpwd=${PWD}

pp5TeVfiles=/rstorage/alice/data/LHC17pq/448/files.txt
pp5TeVMCfiles=/rstorage/alice/data/LHC18b8/449/files.txt
#pp13TeVfiles=/home/software/users/napadula/pyjetty/pyjetty/alihfjets/test_sb_file_list.txt
pp13TeVfiles=/rstorage/napadula/test_hfana/LHC18b_filelist.txt
flist=${pp5TeVMCfiles}

dname=$(date +"%Y-%m-%d-%H-%M")

flistdname=$(dirname ${flist})
#outputdir=$(basename ${flistdname})
outputdir=LHC18b8
outputdir=/rstorage/${USER}/test_hfMC/${outputdir}/${dname}
mkdir -pv ${outputdir}
cd ${outputdir}
pwd

cp -v ${flist} .
split --additional-suffix=.flist -d -l 10 -a 5 ${flist}

job_lists=$(find $PWD -name "*.flist" | sort)

cp -v ${PYJETTY_DIR}/pyjetty/alihfjets/test_hfMC.sh .
cp -v ${PYJETTY_DIR}/pyjetty/alihfjets/test_hfMC.py .

executable=${PWD}/test_hfMC.py

submit_script=${PWD}/submit_all.sh
rm -f ${submit_script}

for fl in ${job_lists}
do
	echo "sbatch --chdir=${PWD} --output=${fl}.output --error=${fl}.error ${PWD}/test_hfMC.sh ${executable} ${fl} ${outputdir}" | tee -a ${submit_script}
	chmod +x ${submit_script}
done

echo "[i] created: ${submit_script}"

cd ${cpwd}
