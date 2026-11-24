import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from matplotlib import animation
from qutip import qfunc
from base64 import b64encode
from IPython.display import HTML
from scipy.ndimage import gaussian_filter

def plot_expect_with_variance(t, op_list, op_title, states):
    """
    Plot the expectation value of an operator (list of operators)
    with an envelope that describes the operators variance.
    """

    fig, axes = plt.subplots(1, len(op_list), figsize=(14, 3))

    for idx, op in enumerate(op_list):
        e_op = qt.expect(op, states)
        v_op = qt.variance(op, states)

        axes[idx].fill_between(
            t, e_op - np.sqrt(v_op), e_op + np.sqrt(v_op), color="green",
            alpha=0.5)
        axes[idx].plot(t, e_op, label="expectation")
        axes[idx].set_xlabel("Time")
        axes[idx].set_title(op_title[idx])

    return fig, axes

def plot_IQnbar_vs_time(t, op, states, name='mode'):

    plot_expect_with_variance(t, [(op+op.dag())/2, (-1j*op+1j*op.dag())/2, op.dag()*op], [name+' I', name+' Q', name+' nbar'], states)
    plt.tight_layout()


def make_Qfunc_movie(t, states, bins=31):
    fig, axes = plt.subplots(1, 1, figsize=(10, 5))

    def display_embedded_video(filename):
        video = open(filename, "rb").read()
        video_encoded = b64encode(video).decode("ascii")
        video_tag = '<video controls alt="test" \
                    src="data:video/x-m4v;base64,{0}">'.format(
            video_encoded
        )
        return HTML(video_tag)

    xvec = np.linspace(-np.sqrt(np.shape(states[0])[0]),-np.sqrt(np.shape(states[0])[0]), bins)

    def update(n):
        print(n)
        axes.cla()
        Q0 = qfunc(states[n], xvec, xvec)
        axes.pcolor(Q0, cmap='magma')
        axes.set_aspect(1)
        return axes.artists
        # return axes[0].artists + axes[1].artists + axes[2].artists

    anim = animation.FuncAnimation(fig, update, frames=np.arange(0, len(states), 1),
                                   blit=True)

    anim.save("saved_animation.mp4", fps=10, writer="ffmpeg")

    plt.close(fig)

    display_embedded_video("saved_animation.mp4")

def Qfunc_fidelity(state0, state1, detector_noise=1):

    lim = int(np.sqrt(np.shape(state0)[0])+1)
    xvec = np.linspace(-lim, lim, 100)

    hist0 = qfunc(state0, xvec, xvec)
    hist1 = qfunc(state1, xvec, xvec)

    filter_sigma = detector_noise / (xvec[1] - xvec[0])

    hist0 = gaussian_filter(hist0, sigma=filter_sigma)
    hist1 = gaussian_filter(hist1, sigma=filter_sigma)

    PDF0 = hist0 / np.sum(hist0)
    PDF1 = hist1 / np.sum(hist1)

    fidelity = 1 - np.sum(np.minimum(PDF0, PDF1)) / 2

    fig, ax = plt.subplots(1, 1, dpi=200, figsize=(3, 3))

    total = (PDF0 + PDF1) / np.max(PDF0 + PDF1)
    img = (PDF1 - PDF0) / (PDF1 + PDF0 + 1e-12)

    G = 0 * img
    R = 0.5 - img / 4
    B = 0.5 + img / 2
    A = total * 1

    rgba = np.stack([R, G, B, A], axis=-1)
    rgba = np.clip(rgba, a_min=0, a_max=1)

    ax.imshow(rgba, extent=[-lim, lim, -lim, lim])

    print('Bayes fidelity: ' + str(fidelity))

    ax.grid()
    ax.set_title('Bayes fidelity: ' + str(fidelity))
    ax.set_xticks(np.linspace(-lim, lim, 5))
    ax.set_yticks(ax.get_xticks())
    plt.show()