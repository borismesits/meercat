import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

def plot_expect_with_variance(N, t, op_list, op_title, states):
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