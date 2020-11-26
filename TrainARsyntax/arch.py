


# from bpnet.losses import multinomial_nll


import keras
import keras.layers as kl
import tensorflow as tf
import tensorflow_probability as tfp
import numpy as np


tasks = ['Oct4', 'Sox2', 'Nanog', 'Klf4']

input = kl.Input(shape=(699, 4))

x = kl.Conv1D(64, kernel_size=25,
    padding='same', activation='relu')(input)

# m1 = keras.models.Model([input], x)

# x = np.random.rand(1,699,4)
#
# m1.predict(x)
#
#
# y = np.random.rand(5,1000,2).astype(np.float32)


for i in range(9):
    conv_x = kl.Conv1D(64, kernel_size=3, padding='same',
        activation='relu', dilation_rate=2**i)(x)
    x = kl.add([conv_x, x])




bottleneck = x


# heads
outputs = []
for task in tasks:
# profile shape head
    px = kl.Reshape((-1, 1, 64))(bottleneck)
    px = kl.Conv2DTranspose(2, kernel_size=(25, 1), padding='same')(px)
    outputs.append(kl.Reshape((-1, 2))(px))
    cx = kl.GlobalAvgPool1D()(bottleneck)
    outputs.append(kl.Dense(2)(cx))




def multinomial_nll(true_counts, logits):
    """Compute the multinomial negative log-likelihood along the sequence (axis=1)
    and sum the values across all each channels
    Args:
      true_counts: observed count values (batch, seqlen, channels)
      logits: predicted logit values (batch, seqlen, channels)
    """
    print(logits)
    # swap axes such that the final axis will be the positional axis
    logits_perm = tf.transpose(logits[0], (0, 2, 1))
    true_counts_perm = tf.transpose(true_counts, (0, 2, 1))
    counts_per_example = tf.reduce_sum(true_counts_perm, axis=-1)
    dist = tfp.distributions.Multinomial(total_count=counts_per_example,
                                                logits=logits_perm)
    # get the sequence length for normalization
    seqlen = float(tf.shape(true_counts)[1])
    return -tf.reduce_sum(dist.log_prob(true_counts_perm)) / seqlen




model = keras.models.Model([input], outputs)
model.compile(keras.optimizers.Adam(lr=0.004),
    loss=[multinomial_nll, 'mse'] * len(tasks),
    loss_weights=[1, 10] * len(tasks))



#############################3


m1 = keras.models.Model([input], x)




x = np.random.rand(5,1000,4)
y = np.random.rand(5,1000,2).astype(np.float32)

logits = model.predict(x)[0]
true_counts = y


logits_perm = tf.transpose(logits, (0, 2, 1))
true_counts_perm = tf.transpose(true_counts, (0, 2, 1))
counts_per_example = tf.reduce_sum(true_counts_perm, axis=-1)
dist = tfp.distributions.Multinomial(total_count=counts_per_example,logits=logits_perm)



log_p = (tf.math.log(self._probs) if self._logits is None else tf.math.log_softmax(self._logits))


dist.log_prob(true_counts_perm)

dtype="float32"
probs = None
tensor_util.convert_nonref_to_tensor(counts_per_example, name='total_count', dtype="float32")

tensor_util.convert_nonref_to_tensor(probs, dtype=dtype, name='probs')

tensor_util.convert_nonref_to_tensor(logits_perm, dtype=dtype, name='logits')


#################333
import tensorflow.compat.v2 as tf

from tensorflow_probability.python.internal import tensor_util
from tensorflow_probability.python.internal import tensorshape_util

@distribution_util.AppendDocstring(_multinomial_sample_note)
def _log_prob(self, counts):
    log_p = (
        tf.math.log(self._probs)
        if self._logits is None else tf.math.log_softmax(self._logits))
    log_unnorm_prob = tf.reduce_sum(
        tf.math.multiply_no_nan(log_p, counts), axis=-1)
    neg_log_normalizer = tfp_math.log_combinations(self.total_count, counts)
    return log_unnorm_prob + neg_log_normalizer

############################


log_p = (
    tf.math.log(probs)
    if logits_perm is None else tf.math.log_softmax(logits_perm))


model.evaluate(x,y)

l = [-50., -43, 0]
d = tfp.distributions.Multinomial(total_count=4., logits=l)







kl.Reshape((-1, 1, 64))(bottleneck)

input => Batch, rows, cols, channels
shape=(None, 1000, 1, 64)
weigth => window1, window2, filter, inpChannel
shape=(25, 1, 2, 64)







t = tf.constant([
                    [
                        [1.0, 1.0, 1.0, 1.0],
                        [2.0, 2.0, 2.0, 2.0],
                        [5.0, 5.0, 5.0, 5.0],
                        [5.0, 5.0, 5.0, 5.0],
                        [5.0, 5.0, 5.0, 5.0],
                        [5.0, 5.0, 5.0, 5.0],
                    ],
                    [
                        [3.0, 3.0, 3.0, 3.0],
                        [4.0, 4.0, 4.0, 4.0],
                        [6.0, 6.0, 6.0, 6.0],
                        [5.0, 5.0, 5.0, 5.0],
                        [5.0, 5.0, 5.0, 5.0],
                        [5.0, 5.0, 5.0, 5.0]
                    ]
                ]
                )


t1 = kl.Conv1D(64, kernel_size=2,    padding='same', activation='relu')(t)

t2 = kl.Reshape((-1, 1, 64))(t1)

t3 = kl.Conv2DTranspose(2, kernel_size=(25, 1), padding='same')(t2)

t4 = kl.Reshape((-1, 2))(t3)



input_shape=(2, 1000, 4)
x = tf.random.normal(input_shape)
y = tf.keras.layers.Conv1D(32, 3, activation='relu',input_shape=input_shape[1:], padding='same')(x)


dx = kl.Conv2DTranspose(2, kernel_size=(25, 1), padding='same')(px)

kl.Reshape((-1, 2))(dx)











np.ones(None, 64)
n = np.array(size)



import numpy as np

tf.matmul(m, np.array([[-0.23118673,  0.22074068],
    [-0.15855451, -0.2774612 ],
    [ 0.24068296, -0.24349943],
    [ 0.17451072,  0.257658  ],
    [-0.00716063, -0.14690058],
    [ 0.25660616,  0.0495317 ],
    [-0.171871  , -0.07320164],
    [ 0.20816362,  0.08312938],
    [ 0.09340459,  0.13587695],
    [-0.24457389, -0.09932217],
    [ 0.08382142, -0.29211196],
    [ 0.07806408,  0.13094017],
    [-0.26552778,  0.09009245],
    [ 0.21853697,  0.08081666],
    [-0.06580997,  0.22820574],
    [ 0.13870049,  0.25198758],
    [-0.23761024,  0.2614932 ],
    [ 0.27793205, -0.229366  ],
    [ 0.13237956,  0.06598672],
    [ 0.00543875,  0.06698093],
    [ 0.24327874, -0.14977501],
    [-0.18502907,  0.26032692],
    [ 0.0660612 ,  0.17036852],
    [-0.11690193,  0.17057252],
    [ 0.0611476 , -0.0957337 ],
    [ 0.13092718,  0.14392787],
    [-0.28531662,  0.12814951],
    [-0.10331975, -0.24110498],
    [ 0.04319856, -0.05830774],
    [ 0.29207438, -0.00319806],
    [-0.13305786, -0.15858412],
    [-0.29131922, -0.05628595],
    [-0.2572455 , -0.23890145],
    [ 0.21256334, -0.03264973],
    [ 0.21419424, -0.21311699],
    [ 0.2027517 , -0.21272385],
    [-0.10803424, -0.10318416],
    [ 0.2184791 ,  0.26885307],
    [-0.15075962, -0.10822998],
    [ 0.11638868, -0.03119725],
    [-0.2975846 , -0.27173662],
    [ 0.18831044, -0.28121486],
    [-0.20391501, -0.14155097],
    [ 0.07232821,  0.10721099],
    [ 0.24559754,  0.1812036 ],
    [-0.04466537,  0.13184255],
    [ 0.12640935,  0.24089271],
    [ 0.07004806,  0.10462123],
    [-0.22069634, -0.10053015],
    [-0.28822792, -0.23324238],
    [-0.22793123, -0.0099692 ],
    [ 0.19485429, -0.00669631],
    [-0.24335471, -0.23557112],
    [-0.25118583,  0.21007335],
    [ 0.21294475,  0.29717165],
    [-0.28491023,  0.1309892 ],
    [ 0.25465202,  0.29358673],
    [ 0.24138671, -0.22947663],
    [ 0.16318926, -0.02272207],
    [-0.09757139,  0.16338357],
    [ 0.08413821,  0.14712372],
    [-0.18241645,  0.25537884],
    [ 0.0052284 ,  0.214858  ],
    [-0.1886361 , -0.24395747]])
    )
