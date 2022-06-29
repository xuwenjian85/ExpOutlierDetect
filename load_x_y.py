file = '/media/eys/xwj/RNAseq/public_normal/myoutput/X_y_20220615-0952.pkl' 
# All 6m sample,gene pairs

print(file)
with open(file, 'rb') as f:
    [ X,y ]=pickle.load(f)
    print(X.shape, y.shape)
