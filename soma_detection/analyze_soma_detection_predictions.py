import numpy as np

y_pred = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/train_test_set/y_pred.npy')
y_pred_2 =  np.load('/data/elowsky/OLST/reconstruction/190327_Emxcre_Reconstruction/soma_detection/model/train_test_set/y_pred.npy')

y_pred = np.concatenate((y_pred,y_pred_2))

y = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/train_test_set/y.npy')
y_2 = np.load('/data/elowsky/OLST/reconstruction/190327_Emxcre_Reconstruction/soma_detection/model/train_test_set/y.npy')

y = np.concatenate((y,y_2))

x_features = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/train_test_set/X_features.npy',allow_pickle=True)


z_unique = np.unique(x_features[:,16])
y_unique = np.unique(x_features[:,17])

fp_nn_out = np.zeros((0,8))
	
##### metrics ######
tp = [1 if a==1 and b == 1 else 0 for a,b in zip(y_pred,y)]
fp = [1 if a==1 and b == 0 else 0 for a,b in zip(y_pred,y)]
fn = [1 if a==0 and b == 1 else 0 for a,b in zip(y_pred,y)]
tn = [1 if a==0 and b == 0 else 0 for a,b in zip(y_pred,y)]
acc_count = [1 if a==b else 0 for a,b in zip(y_pred,y)]

tp_r = []
fp_r = []
fn_r = []
tn_r = []
acc_r = []

if True:
	for i in range(len(tp)):
		#soma = x_features[i]
		#z_f = soma[16]
		if True:
		#if z_f > 0:
			tp_r.append(tp[i])
			fp_r.append(fp[i])
			fn_r.append(fn[i])
			tn_r.append(tn[i])
			acc_r.append(acc_count[i])
		

	tp = tp_r
	fp = fp_r
	fn = fn_r
	tn = tn_r
	acc_count = acc_r

precision = sum(tp) / (sum(tp) + sum(fp))
recall = sum(tp) / (sum(tp) + sum(fn))
f1 = (2*(precision * recall)/(precision + recall))
acc = sum(acc_count) / len(tp)

print('F1:', round(100*f1,1),'%')
print('Precision:',round(100*precision,1),'%')
print('Recall:',round(100*recall,1),'%')
print('Accuracy:',round(100*acc,1),'%')
print()
##################

z_list = [None]*len(z_unique)
y_list = [None]*len(y_unique)


fp_somas = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/1000/fp/fp.npy',allow_pickle=True)
tp_somas = np.load('/data/elowsky/OLST/reconstruction/180206_Emxcre_Reconstruction/soma_detection/model/1000/tp/tp.npy',allow_pickle=True)

for i in range(len(x_features)):

	soma = x_features[i]
	z_f = soma[16]
	y_f = soma[17]

	z_index = np.where(z_unique == z_f)[0][0]
	y_index = np.where(y_unique == y_f)[0][0]

	if z_list[z_index] == None:
		z_list[z_index] = [0,0,0,0]

	if y_list[y_index] == None:
		y_list[y_index] = [0,0,0,0]
	
	z_list[z_index][0] += int(y_pred[i] == 1 and y[i] == 1)
	z_list[z_index][1] += int(y_pred[i] == 1 and y[i] == 0)
	z_list[z_index][2] += int(y_pred[i] == 0 and y[i] == 1)
	z_list[z_index][3] += int(y_pred[i] == 0 and y[i] == 0)

	y_list[y_index][0] += int(y_pred[i] == 1 and y[i] == 1)
	y_list[y_index][1] += int(y_pred[i] == 1 and y[i] == 0)
	y_list[y_index][2] += int(y_pred[i] == 0 and y[i] == 1)
	y_list[y_index][3] += int(y_pred[i] == 0 and y[i] == 0)
	
	if y_pred[i] == 0 and y[i] == 1:
		#soma_number = np.where((soma[0:7] == fp_somas[:,0:7]).all(axis=1))[0][0] + 1
		soma_number = np.where((soma[0:7] == tp_somas[:,0:7]).all(axis=1))[0][0] + 1
		fp_nn_out = np.vstack((fp_nn_out,np.hstack((soma[0:7],soma_number))))

precision_z = [-1 if (x[0] + x[1]) == 0 else int(100*x[0] / (x[0] + x[1])) for x in z_list]
recall_z = [-1 if (x[0] + x[2]) == 0 else int(100*x[0] / (x[0] + x[2])) for x in z_list]

precision_y = [-1 if (x[0] + x[1]) == 0 else int(100*x[0] / (x[0] + x[1])) for x in y_list]
recall_y = [-1 if (x[0] + x[2]) == 0 else int(100*x[0] / (x[0] + x[2])) for x in y_list]

print(z_unique)
print(precision_z)
print(recall_z)
print(z_list)
fp_nn_out = fp_nn_out[fp_nn_out[:,3].argsort()]
#np.savetxt('/home/elowsky/Desktop/fn_nn_out.csv',fp_nn_out,delimiter=",",fmt="%s")






