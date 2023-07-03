
#How to end test
def end_test(figures,suffix=''):
    print(f"Saving image for {suffix}...")
    for ifig,fig in enumerate(figures):
        fig.savefig(f"/tmp/test{suffix}-{ifig:03d}.png")
    plt.close("all")
    return
