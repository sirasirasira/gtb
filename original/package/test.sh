<<commentout

export DYLD_LIBRARY_PATH=/Users/takigawa/Desktop/graphlearn/boost_local/lib/:$DYLD_LIBRARY_PATH
echo glearn
time ./glearn/glearn -1 0.03 ./sample_data/cpdb.train --outfile=./out/model1.out > /dev/null

echo glearn2
time ./glearn2/glearn -1 0.03 ./sample_data/cpdb.train --outfile=./out/model2.out > /dev/null

echo training error 1
./eval/evaluator ./out/model1.out ./sample_data/cpdb.train | tail -n 1

echo test error 1
./eval/evaluator ./out/model1.out ./sample_data/cpdb.test | tail -n 1

echo training error 2
./eval/evaluator ./out/model2.out ./sample_data/cpdb.train | tail -n 1

echo test error 2
./eval/evaluator ./out/model2.out ./sample_data/cpdb.test | tail -n 1

echo generate subgraph existence
python ./py/extract_pat.py ./out/model1.out > ./out/test_feat.out
./finder/finder --tran_by_feat=./out/cpdb_train.mat --feat_by_tran=/dev/null ./out/test_feat.out ./sample_data/cpdb.train > /dev/null
./finder/finder --tran_by_feat=./out/cpdb_test.mat --feat_by_tran=/dev/null ./out/test_feat.out ./sample_data/cpdb.test > /dev/null

echo test machine learning
s=`cat out/test_feat.out | wc -l`
python py/machine_learning_sample.py ./out/cpdb_train.mat ./out/cpdb_test.mat $s

commentout

echo generate frequent subgraphs
./gspan/gspan -x 5 ./sample_data/cpdb.train > ./out/freqpat_maxpat5.out
python ./py/extract_pat.py ./out/freqpat_maxpat5.out > ./out/gf_maxpat5_feat.out
./finder/finder --tran_by_feat=./out/cpdb_train_gf.mat --feat_by_tran=/dev/null ./out/gf_maxpat5_feat.out ./sample_data/cpdb.train > /dev/null
./finder/finder --tran_by_feat=./out/cpdb_test_gf.mat --feat_by_tran=/dev/null ./out/gf_maxpat5_feat.out ./sample_data/cpdb.test > /dev/null

echo test machine learning
s=`cat out/gf_maxpat5_feat.out | wc -l`
#python py/machine_learning_sample.py ./out/cpdb_train_gf.mat ./out/cpdb_test_gf.mat $s
python py/gradient_boosting.py ./out/cpdb_train_gf.mat ./out/cpdb_test_gf.mat $s

