
sessiondir=/data/ppp/00.PredictionOnly01
query=/data/ppp/00.PredictionOnly01/query
modelName=/data/ppp/models/Parageobacillus_thermoglucosidasius.cnn_lstm_71/cnn_lstm.h5


python3 /data/ppp/ppp_Prediction_Only_Anne.py -sessiondir $sessiondir -query $query -modelName $modelName
