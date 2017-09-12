for i in `seq 0 9`;
do
  ./corels -c 2 -v 10 -n 100000000 -p 1 -r 0.005 ../data/CrossValidation/compas_"$i"_train.out ../data/CrossValidation/compas_"$i"_train.label ../data/CrossValidation/compas_"$i"_train.minor
done

