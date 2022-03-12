import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat,savemat
import os
# Tools for saving images
import scipy as sp
import scipy.misc
import imageio
import PIL
from PIL import Image
Image.MAX_IMAGE_PIXELS
Image.MAX_IMAGE_PIXELS=1e10 # forget attack
import tifffile as tiff

###############################################################
# Helper Functions

def inhomogeneity_correction_2d(I,dI,Nbins,a,p,sigmaR,niter,ep,expand_factor=0.5,return_figure=False):
    '''
    Correct inhomogeneity in a 2D image with a smooth multiplicative bias field.
    
    Parameters
    ----------
    I : numpy array
        Input image to be corrected
    dI : 2 element iterable
        Pixel size along each dimension
    Nbins : int
        Number of bins for kernel density estimate computation. Should enough that 
        suming over bins is a good approximation to an integral.
    a : float
        Spatial scale of smoothing. Highpass operator: (id - a^2 Delta)^p
    p : float
        Power of smoothing. Highpass operator: (id - a^2 Delta)^p
    sigmaR : float
        Regularization strength of the form 1/2/sigmaR^2
    niter : int
        Number of iterations of gradient descent
    ep : float
        Step size for gradient descent
    expand_factor : float
        Expand historam away from min and max value of image by this factor. Default 0.5.
    return_figure : bool
        If true return a figure in output arguments.
    
    Returns
    -------
    J : numpy array
        Corrected image
    u : numpy array
        Bias field for applying correction: J = I * exp(u)
    f : matplotlib figure handle
        Figure showing results if return_figure is True.
    '''
    
    # setup variables
    minI = np.min(I)
    maxI = np.max(I)
    diffI = maxI - minI

    mint = minI - diffI*expand_factor
    maxt = maxI + diffI*expand_factor

    bins = np.linspace(mint,maxt,Nbins)
    delta = bins[1]-bins[0]
    sigma = delta*0.5

    # initialize
    Ipad = np.concatenate([np.concatenate([I,I[:,::-1]],1), np.concatenate([I[::-1],I[::-1,::-1]],1)],0)
    u = np.zeros_like(Ipad)

    # smoothing
    nI = np.array(Ipad.shape)
    xI = [np.arange(n)*d for n,d in zip(nI,dI)]
    fI = [np.arange(n)/n/d for n,d in zip(nI,dI)]
    FI = np.stack(np.meshgrid(*fI,indexing='ij'))

    L = (1.0 - 2.0*a**2*(  (np.cos(2.0*np.pi*FI[0]) - 1.0)/dI[0]**2 + (np.cos(2.0*np.pi*FI[1]) - 1.0)/dI[1]**2  ))**p
    LL = L**2
    K = 1.0/LL 


    
    # run the algorithm
    f,ax = plt.subplots(2,3)
    Esave = []
    regsave = []
    Hsave = []

    for it in range(niter):
        # the corrected image
        J = Ipad * np.exp(u)

        # get the probability density and gradient    
        pdf = np.zeros_like(bins)
        grad = np.zeros_like(u)
        for i in range(len(bins)):
            k =  np.exp(-0.5*( bins[i] - J )**2/sigma**2)/np.sqrt(2.0*np.pi*sigma**2) 
            pdf[i] = np.sum( k )/I.size
            grad += (np.log(pdf[i]) + 1)*delta/J.size *  k * (J - bins[i])/sigma**2 * J 

        # entropy
        H = -np.sum(pdf*np.log(pdf))*delta
        # reg
        reg = np.sum(np.abs(np.fft.fftn(u))**2*LL)/u.size*np.prod(dI)/2/sigmaR**2    
        # cost
        E = H + reg
        # save
        Hsave.append(H)
        regsave.append(reg)
        Esave.append(E)

        # smooth it
        grad = np.fft.ifftn(np.fft.fftn(grad)*K).real
        # add reg 
        grad += u/sigmaR**2*np.prod(dI)

        # update
        u -= grad*ep

        # visualize
        '''
        # draw the original image
        ax[0,0].cla()
        ax[0,0].imshow(I,cmap='gray')
        ax[0,0].set_title('Original')
        # draw the corrected image
        ax[1,0].cla()
        ax[1,0].imshow(J[:I.shape[0],:I.shape[1]],cmap='gray')
        ax[1,0].set_title('Corrected')
        # draw the cost
        ax[0,2].cla()
        ax[0,2].plot(Esave)
        ax[0,2].set_title('objective')
        ax[1,1].cla()
        ax[1,1].plot(regsave)
        ax[1,1].set_title('Regularization')
        ax[1,2].cla()
        ax[1,2].plot(Hsave)
        ax[1,2].set_title('Entropy')
        # correction
        ax[0,1].cla()
        ax[0,1].imshow(u[:I.shape[0],:I.shape[1]])
        ax[0,1].set_title('Bias field')

        # render
        f.canvas.draw()
        '''
    J = J[:I.shape[0],:I.shape[1]]
    u = u[:I.shape[0],:I.shape[1]]
    output = (J,u)
    if return_figure:
        output = (J,u,f)
    return output

#################################################################
# Wrapper Functions 

def runBias(brainNum,blockNum,mriFolder,mriSuffix,fname):
    data = loadmat(fname)
    
    # inputs
    # image
    I = data['mris'][0][2]
    # pixel size
    dI = np.array([1.0,1.0]) #

    # number of bins
    Nbins = 20
    # how much to expand range relative to max-min
    expand_factor = 0.5
    # smoothness, highpass penalty L = (1 - a^2 Delta)^p
    # penality is 1/2/sigmaR^2 |Lu|^2_{L2}
    a = 10.0
    p = 2.0
    sigmaR = 5e2

    niter = 200
    ep = 2e4


    I = data['mris'][0][0]
    data['corrected'] = [[]]
    data['bias'] = [[]]

    names = ['L1','L5','L10','L15','L20','L25','L30','L35','L40','L45','L50','L55','L60','L65','L70','L75']
    if (brainNum == 2):
        names = ['L5','L6','L7','L8','L9','L10','L11','L12','L13','L14','L15','L16','L17','L18','L19']
    elif (brainNum == 5 and blockNum == 3):
        names = ['L3','L7','L12','L17','L22','L27','L32','L37']
    num = len(data['mris'][0])
    for i in range(num):
        I = data['mris'][0][i]
        J,u,f = inhomogeneity_correction_2d(I,dI,Nbins,
                                            a,p,sigmaR,
                                            niter,ep,
                                            expand_factor=0.5,return_figure=True)
        data['corrected'][0].append(J)
        data['bias'][0].append(u)
        f.savefig(f'output_{i:03d}.jpg')
        f.canvas.draw()
        fig,ax = plt.subplots()
        ax.imshow(J,cmap='gray')
        fig.canvas.draw()
        fig.savefig('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/histology/alignment/block' + str(blockNum) + '/' + mriFolder + '/' + names[i]+ mriSuffix + '.png')
        print(J.dtype)
        print("max and min of J " + str(np.min(J)) + ", " + str(np.max(J)))
        Jrange = np.max(J) - np.min(J)
        Jnorm = (J - np.min(J))/Jrange
        J = Jnorm.astype('float32')
        #J = J.astype('uint8')
        #tiff.imsave('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain5/histology/alignment/block2/entropy/' + names[i]+'.tif',J)
        im = Image.fromarray(J)
        im.save('/cis/home/kstouff4/Documents/datasets/exvivohuman_11T/more_blocks/Brain' + str(brainNum) + '/histology/alignment/block' + str(blockNum) + '/' + mriFolder + '/' + names[i]+ mriSuffix+ '.tif')
        J_ = (J - np.mean(J))/np.std(J)*np.std(I) + np.mean(I)
       
    return
