#!/bin/bash

cpwd=${PWD}

pp5TeVfiles=/rstorage/alice/data/LHC17pq/448/files.txt
pp5TeVMCfiles=/rstorage/alice/data/LHC18b8/449/files.txt
pp13TeVfiles=/rstorage/napadula/test_hfana/LHC18l_filelist.txt
flist=${pp13TeVfiles}

dname=$(date +"%Y-%m-%d-%H-%M")

flistdname=$(dirname ${flist})
#outputdir=$(basename ${flistdname})
outputdir=LHC18l
outputdir=/rstorage/${USER}/test_event/${outputdir}/${dname}
mkdir -pv ${outputdir}
cd ${outputdir}
pwd

cp -v ${flist} .
split --additional-suffix=.flist -d -l 15 -a 5 ${flist}

job_lists=$(find $PWD -name "*.flist" | sort)

cp -v ${PYJETTY_DIR}/pyjetty/alihfjets/test_event.sh .
cp -v ${PYJETTY_DIR}/pyjetty/alihfjets/test_event.py .

executable=${PWD}/test_event.py

submit_script=${PWD}/submit_all.sh
rm -f ${submit_script}

for fl in ${job_lists}
do
	echo "sbatch --chdir=${PWD} --output=${fl}.output --error=${fl}.error ${PWD}/test_event.sh ${executable} ${fl} ${outputdir}" | tee -a ${submit_script}
	chmod +x ${submit_script}
done

echo "[i] created: ${submit_script}"

cd ${cpwd}
