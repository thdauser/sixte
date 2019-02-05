#! /bin/bash

script_name=run_test
log_file=run_e2e_test.log


if [ $# -eq 0 ]
then
    list_dirs=( $(ls -d *.e2e) )
elif [ $# -eq 1 ]
then     
     list_dirs=$1
else    
    echo "Usage: $0 [dir] "
    exit 1
fi

n_tests=0
n_succ=0
n_fail=0

function run_e2e_test_single {
    name_dir=$1
    n_total=$2
    cd $name_dir

    n_tests=$(($n_tests+1))

    echo -ne " - running e2e test in $name_dir ($n_tests/$n_total) \t " |tee $log_file

    if  ./$script_name >$name_dir.log 2>$name_dir.err ; then
	echo " SUCCESS "
	n_succ=$(($n_succ+1))
    else
	echo " FAILURE "
	n_fail=$(($n_fail+1))
    fi
}

echo " "
echo " *****  Performing E2E Testing *****"  |tee $log_file
echo " "  |tee $log_file

getcwd=`pwd`

n_total=${#list_dirs[@]}
for name_dir in ${list_dirs[@]}; do
        
    if [ -e $name_dir/$script_name ]
    then
	run_e2e_test_single $name_dir $n_total
    else
	echo " ** warning: no script called $script_name in $name_dir -> skipping test "
    fi

    # need to get back to base dir
    cd $getcwd
done

echo "  "  |tee $log_file
echo " *****  SUMMARY of E2E Testing *****"  |tee $log_file
echo "     TOTAL   : $n_tests"  |tee $log_file
echo "     SUCESS  : $n_succ" |tee $log_file
echo "     FAILED  : $n_fail" |tee $log_file 
echo " ***********************************"  |tee $log_file

if (($n_fail>0)) ; then
echo "  End-to-End test FAILED "  |tee $log_file
exit 1
else
echo "  End-to-End test SUCCESSFUL "  |tee $log_file
fi

echo "  "  |tee $log_file
