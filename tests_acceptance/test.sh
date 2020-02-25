set -euo pipefail 
echo "Starting tests"
../src/corels -k 213 -n 100000002 -t 1 -v 10 -c 2 -p 1 -r 0.01 ../data/CrossValidation/compas_1_train.out ../data/CrossValidation/compas_1_train.label ../data/CrossValidation/compas_1_train.minor
../src/corels -k 9031 -n 100000002 -t 4 -v 10 -c 2 -p 1 -r 0.01 ../data/CrossValidation/compas_1_train.out ../data/CrossValidation/compas_1_train.label ../data/CrossValidation/compas_1_train.minor
../src/corels -k 82 -n 100000002 -t 8 -v 10 -c 2 -p 1 -r 0.01 ../data/CrossValidation/compas_1_train.out ../data/CrossValidation/compas_1_train.label ../data/CrossValidation/compas_1_train.minor
echo "Compas done"
../src/corels -k 1902 -n 100000002 -v 1 -c 2 -p 1 -t 1 -r 0.01 ../data/haberman.out ../data/haberman.label ../data/haberman.minor
../src/corels -k 200 -n 100000002 -v 1 -c 2 -p 1 -t 4 -r 0.01 ../data/haberman.out ../data/haberman.label ../data/haberman.minor
../src/corels -k 783 -n 100000002 -v 1 -c 2 -p 1 -t 8 -r 0.01 ../data/haberman.out ../data/haberman.label ../data/haberman.minor
echo "Haberman done"
../src/corels -k 651 -t 1 -c 2 -n 10000000 -p 1 -v 10 -r 0.01 ../data/CrossValidation/propublica_ours_1_train.out ../data/CrossValidation/propublica_ours_1_train.label ../data/CrossValidation/propublica_ours_1_train.minor
../src/corels -k 823 -t 4 -c 2 -n 10000000 -p 1 -v 10 -r 0.01 ../data/CrossValidation/propublica_ours_1_train.out ../data/CrossValidation/propublica_ours_1_train.label ../data/CrossValidation/propublica_ours_1_train.minor
../src/corels -k 22 -t 8 -c 2 -n 10000000 -p 1 -v 10 -r 0.01 ../data/CrossValidation/propublica_ours_1_train.out ../data/CrossValidation/propublica_ours_1_train.label ../data/CrossValidation/propublica_ours_1_train.minor
echo "Propublica done"
../src/corels -k 2672 -t 1 -c 2 -p 1 -r 0.01 -n 2000000000 ../data/1-120-5000-0.out ../data/1-120-5000-0.label ../data/1-120-5000-0.minor 
../src/corels -k 262 -t 4 -c 2 -p 1 -r 0.01 -n 2000000000 ../data/1-120-5000-0.out ../data/1-120-5000-0.label ../data/1-120-5000-0.minor 
../src/corels -k 72 -t 8 -c 2 -p 1 -r 0.01 -n 2000000000 ../data/1-120-5000-0.out ../data/1-120-5000-0.label ../data/1-120-5000-0.minor 
echo "Census done"
echo "All done"