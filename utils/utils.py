# additional functions
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
import numpy as np
import torch

def plot_root_intervals(w, b, interval_of_search=(-20, 60), eps=0.01):
    """
    This function plots intervals of ages where monomials do not exceed probability interval. 
    eps = 0.01 - small constant to not touch to the interval boundaries where p=0 or p=1
    """
    x_intervals = np.array([[-bi/wi if -bi/wi < (1-bi)/wi else (1-bi)/wi for bi, wi in zip(b, w)], 
                            [(1-bi)/wi if -bi/wi < (1-bi)/wi else -bi/wi for bi, wi in zip(b, w)],]).T
    sorted_intervals = x_intervals[np.argsort(x_intervals[:, 0])]
    interval_of_search_corrected = (x_intervals[:, 0].max() + eps, x_intervals[:, 1].min() - eps)
    print('Corrected interval:', interval_of_search_corrected)

    x = np.linspace(interval_of_search[0], interval_of_search[1], 100)
    share_of_p0 = []
    share_of_p1 = []
    for xi in x:
        p = w * xi + b
        share_of_p0.append((p < 0).sum() / p.shape[0])
        share_of_p1.append((p > 1).sum() / p.shape[0])

    fig, axes = plt.subplots(1, 4, figsize=(20, 4))
    for i, (xmin, xmax) in enumerate(sorted_intervals):
        axes[0].hlines(i, xmin, xmax)
    axes[0].axvline(interval_of_search[0], lw=1, ls='--', color='red')
    axes[0].axvline(interval_of_search[1], lw=1, ls='--', color='red')
    axes[0].set_title('Intervals where 0<p<1')
    axes[0].set_xlabel('Age')
    axes[0].set_ylabel('CpG id');

    for i, (xmin, xmax) in enumerate(sorted_intervals):
        axes[1].hlines(i, xmin, xmax)
    axes[1].axvline(interval_of_search[0], lw=2, ls='--', color='red')
    axes[1].axvline(interval_of_search[1], lw=2, ls='--', color='red')
    axes[1].axvline(interval_of_search_corrected[0], lw=2, ls='--', color='orange')
    axes[1].axvline(interval_of_search_corrected[1], lw=2, ls='--', color='orange')
    axes[1].set_title('Intervals where 0<p<1 (detailed)')
    axes[1].set_xlim([interval_of_search[0] - 5, interval_of_search[1] + 5])
    axes[1].set_xlabel('Age')
    axes[1].set_ylabel('CpG id');

    axes[2].plot(x, share_of_p0)
    axes[3].plot(x, share_of_p1)
    axes[2].axvline(interval_of_search_corrected[0], lw=2, ls='--', color='orange')
    axes[2].axvline(interval_of_search_corrected[1], lw=2, ls='--', color='orange')
    axes[3].axvline(interval_of_search_corrected[0], lw=2, ls='--', color='orange')
    axes[3].axvline(interval_of_search_corrected[1], lw=2, ls='--', color='orange')
    axes[2].set_title('Share of p~0')
    axes[3].set_title('Share of p~1')
    axes[2].set_xlabel('Age')
    axes[3].set_xlabel('Age');
    plt.tight_layout()
    plt.show()


def get_max(w, b, domain=None, verbose=0, proper_interval=True):
    """
    Get maxium of Likelihood by checking each concave part of polynimial
    likelihood separately.
    
    Params:
    w - weights of CpGs from reference data
    b - biases of CpGs from reference data
    domain - interval of search
    verbose - if print the information about whether the solution was on the boundary of the interval
    proper_interval - if domain includes ages where p < 0 or p > 1, then do not use log-likelihood
    and use direct product instead

    Return: optimal age
    """
    roots = -b / w
    ord = 'even' if (len(roots) % 2) == 0 else 'odd'
    sign = np.prod(np.sign(w))
    roots = np.sort(roots)
    domain = [roots[0], roots[-1]] if domain is None else domain

    if proper_interval:
        ff = lambda x: np.sum([np.log(wi*x + bi) for wi, bi in zip(w, b)])
    else:
        ff = lambda x: np.prod([(wi*x + bi) for wi, bi in zip(w, b)])
    df = lambda x: np.sum([wi/(wi*x + bi) for wi, bi in zip(w, b)])

    if (ord == 'even') and (sign > 0): #first root left from min; last root right from min
        pairs = roots[1:-1].reshape(-1, 2).tolist()
    elif (ord == 'even') and (sign < 0): #first root left from max; last root right from max
        pairs = roots.reshape(-1, 2).tolist()
    elif (ord == 'odd') and (sign > 0): #first root left from max; last root right from min
        pairs = roots[:-1].reshape(-1, 2).tolist()
    elif (ord == 'odd') and (sign < 0): #first root left from min; last root right from max
        pairs = roots[1:].reshape(-1, 2).tolist()
        
    maxroots = []
    for pair in pairs:
        if pd.Interval(*domain).overlaps(pd.Interval(*pair)):
            root = fsolve(df, np.mean(pair))[0]
        else:
            continue
        if (domain[0] <= root <= domain[1]): #maximum must be within the interval
            maxroots.append(root)
    
    maxroots = [domain[0]] + maxroots + [domain[1]] # check interval boundaries
    solution = maxroots[np.argmax([ff(m) for m in maxroots])]
    if verbose > 0:
        if (solution == maxroots[0]):
            print('Solution on the start of the interval!')
        elif solution == maxroots[-1]:
            print('Solution on the end of the interval!')
        else:
            pass
    return solution


class BetaModel(torch.nn.Module):
    def __init__(self, w, b, t0):
        super(BetaModel, self).__init__()
        self.w = w
        self.b = b
        self.t = torch.nn.Parameter(t0, requires_grad=True)
        self.logsigm = torch.nn.LogSigmoid() 

    def forward(self, c):
        g = self.w * self.t + self.b
        return (c * self.logsigm(g) + (1 - c) * self.logsigm(1 - g))
