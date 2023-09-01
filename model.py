

import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from torch.autograd import Variable
import torch.nn.functional as F
import os
from google.colab import drive
from sklearn.metrics import r2_score

use_cuda = torch.cuda.is_available()
device = torch.device("cuda:0" if use_cuda else "cpu")
torch.backends.cudnn.benchmark = True

def seed_torch(seed):
    import random
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if device.type == 'cuda':
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed) # if you are using multi-GPU.
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

seed_torch(1)

drive.mount('/content/drive')

targetData = pd.read_csv('/content/drive/MyDrive/TestData/target_log2fc.csv')

countsData = pd.read_csv('/content/drive/MyDrive/TestData/counts_log2fc.csv')

targetData = targetData.drop('Unnamed: 0', axis=1)
targetData.head()

countsData = countsData.drop('Unnamed: 0', axis=1)
countsData.head()

countsDataframe = countsData
countsDataframe.head()

tempdf = pd.concat([countsData, targetData], axis=1)
tempdf = tempdf.dropna()

from scipy import stats
tempdf_filtered = tempdf[(np.abs(stats.zscore(tempdf['V1'].to_frame())) < 3).all(axis=1)]

targetData = tempdf_filtered['V1'].to_frame()
countsData = tempdf_filtered.drop(columns='V1')

targetData[targetData < 0] = 0
targetData[targetData > 0] = 1

countsData = countsData.to_numpy()
targetData = targetData.to_numpy()

from sklearn import preprocessing
scalerCounts = preprocessing.StandardScaler().fit(countsData)
scalerTarget = preprocessing.StandardScaler().fit(targetData)

countsData = scalerCounts.transform(countsData)

countsData = countsData.astype(np.float32)
targetData = targetData.astype(np.float32)

counts_tensor = torch.from_numpy(countsData)
target_tensor = torch.from_numpy(targetData)
target_tensor = target_tensor.flatten()

'''
from torch.nn.functional import normalize

counts_tensor = normalize(counts_tensor, p=1.0, dim=1)
counts_tensor
'''

target_tensor.shape

torch.min(target_tensor), torch.max(target_tensor)

target_tensor.shape

#projection.shape

class CustomDataset(torch.utils.data.Dataset):
  def __init__(self,cellstart,cellend):
    self.cellidx=[cellstart+i for i in range(cellend-cellstart+1)]

  def __len__(self):
    return len(self.cellidx)

  def __getitem__(self,index):
    cidx=self.cellidx[index]
    outputrna=counts_tensor[cidx,:]
    outputtarget=target_tensor[cidx]
    return outputrna,outputtarget

import random
def seed_worker(worker_id):
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_seed)
    random.seed(worker_seed)

g = torch.Generator()
g.manual_seed(0)

num_workers=0

loss = nn.CrossEntropyLoss()
U, S, V = torch.pca_lowrank(counts_tensor[0:1060,], q=40)
Vs = V.detach()
V=V.to(device)
V.shape
#approx = torch.matmul(projection, Vs.T)
#mse=loss(approx,counts_tensor)

from torch.utils.data import DataLoader
traindata=CustomDataset(0,1060)
testdata=CustomDataset(1061,1326)

train_dataloader = DataLoader(traindata, batch_size=32, shuffle=True, num_workers=num_workers, worker_init_fn=seed_worker, generator=g)
test_dataloader = DataLoader(testdata, batch_size=32, num_workers=num_workers, worker_init_fn=seed_worker, generator=g)

from torch.nn.modules.linear import Linear
from torch.nn.modules.activation import ReLU
from torch.nn.modules.batchnorm import BatchNorm1d
#model
class Net(nn.Module):
  def __init__(self):
    super(Net, self).__init__()
    self.fc_layer = nn.Sequential(

        nn.Linear(40, 500),
        nn.Dropout(),
        nn.BatchNorm1d(500),
        nn.ReLU(),
        nn.Linear(500, 250),
        nn.BatchNorm1d(250),
        nn.ReLU(),

        nn.Linear(250, 125),
        nn.Dropout(),
        nn.BatchNorm1d(125),
        nn.ReLU(),
        nn.Linear(125, 60),
        nn.BatchNorm1d(60),
        nn.ReLU(),

        nn.Linear(60, 30),
        nn.Dropout(),
        nn.BatchNorm1d(30),
        nn.ReLU(),
        nn.Linear(30, 10),
        nn.BatchNorm1d(10),
        nn.ReLU(),


        nn.Linear(10, 2)

    )
  def forward(self, x):

    x = self.fc_layer(x)

    return x

model = Net().to(device)

num_epochs=100
learning_rate = 0.01

criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=1e-3)

test_loss = []
train_loss = []
test_counter = []
train_counter = []
predicted_values = []
actual_values = []

def cust_loss(output, copy):
  loss = (output - copy).sum()
  return loss

def train(epoch):
  model.train()
  total_step = len(train_dataloader)
  total_train_loss = 0
  for idx, (outputrna, outputtarget) in enumerate(train_dataloader):
        # Move tensors to the configured device
        outputtarget = outputtarget.type(torch.LongTensor)
        outputrna = outputrna.to(device)
        outputtarget = outputtarget.to(device)
        projection = torch.matmul(outputrna, V)
        projection = projection.to(device)

        # Forward pass
        outputs = model(projection)
        loss = criterion(outputs, outputtarget)


        # Backpropagation and optimization
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        total_train_loss += loss.item()
        if (idx+1) % 1 == 0:
            print ('Epoch [{}/{}], Step [{}/{}], Loss: {:.4f}'
                   .format(epoch, num_epochs, idx, total_step, loss.item()))
  total_train_loss = total_train_loss / (idx + 1)
  train_loss.append(total_train_loss)

def test():
  model.eval()
  total_test_loss = 0
  with torch.no_grad():
    correct = 0
    total = 0
    for idx, (outputrna, outputtarget) in enumerate(test_dataloader):
        outputtarget = outputtarget.type(torch.LongTensor)
        outputrna = outputrna.to(device)
        outputtarget = outputtarget.to(device)
        projection = torch.matmul(outputrna, V)
        projection = projection.to(device)
        outputs = model(projection)
        _, predicted = torch.max(outputs.data, 1)
        total += outputtarget.size(0)
        correct += (predicted == outputtarget).sum().item()
        loss = criterion(outputs, outputtarget)
        total_test_loss += loss.item()

    print('Accuracy of the network: {} %'.format(100 * correct / total))
    total_test_loss = total_test_loss / (idx + 1)
    test_loss.append(total_test_loss)

gradient_list = []

def gradient():
  model.eval()
  correct = 0
  total = 0
  for idx, (outputrna, outputtarget) in enumerate(test_dataloader):
      outputtarget = outputtarget.type(torch.LongTensor)
      outputrna = outputrna.to(device)
      outputtarget = outputtarget.to(device)

      outrna=outputrna.detach().requires_grad_(True)
      projection = torch.matmul(outrna, V)
      projection = projection.to(device)
      outputs = model(projection)
      _, predicted = torch.max(outputs.data, 1)
      total += outputtarget.size(0)
      correct += (predicted == outputtarget).sum().item()
      loss = criterion(outputs, outputtarget)

      #simple gradient
      output_copy = outputs.detach().clone()
      loss1 = cust_loss(outputs, output_copy)
      loss1.backward()
      gradient_list.append(outrna.grad)

test()
for epoch in range(1, num_epochs + 1):
  train(epoch)
  test()
gradient()

grd = torch.cat(gradient_list, dim=0).to(device)
grd.shape

relu = nn.ReLU()

cMatrix = counts_tensor[1061:1327, :].to(device)

grad_score = torch.mul(cMatrix, relu(grd))

grad_score_simp = torch.mul(cMatrix, relu(grd*(-1)))

grad_scorenp = grad_score.to('cpu').detach().numpy().copy()
grad_score_simpnp = grad_score_simp.to('cpu').detach().numpy().copy()
grad_scorenp

grad_scoredf = pd.DataFrame(grad_scorenp)
grad_scoredf
grad_simp_scoredf = pd.DataFrame(grad_score_simpnp)

grad_scoredf.columns = countsDataframe[1061:1327].columns
grad_simp_scoredf.columns = countsDataframe[1061:1327].columns
grad_scoredf

import statistics
grad_scored_means = pd.DataFrame(index = ['Means'], columns = grad_scoredf.columns)
for i in range(len(grad_scoredf.columns)):
  mean = statistics.mean(grad_scoredf.iloc()[:,i].sort_values(ascending = False).head(50))
  grad_scored_means.iloc()[:,i] = mean

grad_scored_means

import statistics
grad_scored_means_simp = pd.DataFrame(index = ['Means'], columns = grad_simp_scoredf.columns)
for i in range(len(grad_simp_scoredf.columns)):
  mean = statistics.mean(grad_simp_scoredf.iloc()[:,i].sort_values(ascending = False).head(50))
  grad_scored_means_simp.iloc()[:,i] = mean

grad_scored_means_simp

grad_scored_means_bot_t = grad_scored_means_simp.T
grad_bot150 = grad_scored_means_bot_t.sort_values(by = "Means", ascending = False).head(150)
grad_bot150

grad_scored_means_t = grad_scored_means.T
grad_top150 = grad_scored_means_t.sort_values(by = "Means", ascending = False).head(150)
grad_top150

top50genes = grad_top150.index.tolist()
top50genes
top50genes_str = ','.join(top50genes)
top50genes_str

import statistics
grad_scored_std = pd.DataFrame(index = ['STD'], columns = grad_scoredf.columns)
for i in range(len(grad_scoredf.columns)):
  stdev = statistics.stdev(grad_scoredf.iloc()[:,i])
  grad_scored_std.iloc()[:,i] = stdev

grad_scored_std

grad_scored_std_t = grad_scored_std.T
gradstd_top50 = grad_scored_std_t.sort_values(by = "STD", ascending = False).head(50)
gradstd_top50

testL, = plt.plot(test_loss, label='Test Loss')
trainL, = plt.plot(train_loss, label='Train Loss')
leg = plt.legend(loc='upper center')
plt.ylim([0, 1])
plt.show()
