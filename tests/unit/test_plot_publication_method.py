import matplotlib.pyplot as plt

from benchmark_tools.plot_publication_method import EDGES, NODES, plot


def test_workflow_order_and_separate_outputs():
    names = [node[0] for node in NODES]
    assert len(names) == len(set(names)) == 10
    assert names.index("profiles") < names.index("groups") < names.index("candidates")
    assert names.index("reconcile") < names.index("constraints") < names.index("outputs")
    assert EDGES == tuple(zip(names[:-1], names[1:]))
    figure = plot()
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    boxes = [text.get_window_extent(renderer) for text in figure.axes[0].texts if text.get_text()]
    frame = figure.bbox
    for box in boxes:
        assert frame.contains(box.x0, box.y0) and frame.contains(box.x1, box.y1)
    for i, box in enumerate(boxes):
        assert not any(box.overlaps(other) for other in boxes[i + 1:])
    assert len(figure.axes[0].patches) == 11
    plt.close(figure)
