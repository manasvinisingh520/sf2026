from tensorflow.core.util import event_pb2
from tensorflow.python.lib.io import tf_record
from pathlib import Path

def read_tensorboard_events(logdir):
    """Read all events from a logdir, yield (step, tag, value) tuples."""
    for f in Path(logdir).rglob("events.out.tfevents*"):
        for record in tf_record.tf_record_iterator(str(f)):
            event = event_pb2.Event.FromString(record)
            for value in event.summary.value:
                if value.HasField("simple_value"):
                    yield event.step, value.tag, value.simple_value

def check_overfitting(logdir, tag="Loss/Test", margin=0.01, min_vals=20):
    """Return 1 if overfitting flagged, 0 otherwise. Needs >= min_vals Loss/Test values."""
    loss_values = []
    for step, t, value in read_tensorboard_events(logdir):
        if t == tag:
            loss_values.append(value)

    if len(loss_values) < 200:
        print (f"Not enough loss values to check for overfitting in {logdir}: {len(loss_values)}")

    # Run running average of window of size 20
    running_average = []
    for i in range(len(loss_values)):
        if i < 20:
            running_average.append(loss_values[i])
        else:
            running_average.append(sum(loss_values[i-20:i]) / 20)

    loss_values = running_average

    last_n = loss_values[-min_vals:]
    min_at_end = min(last_n)
    min_overall = min(loss_values)
    min_overall_index = loss_values.index(min_overall)
    if min_at_end > (min_overall + margin) and min_overall_index < 100:
        return 1
    else:
        return 0

# Fuction to count overfitts across multiple testfolds
def count_overfitts(baseDirPath = None):
    """Count the number of overfitts in a logdir."""
    overfitts = 0
    subdirs = [f"{baseDirPath}_testfold{i}" for i in range(1, 9)]
    for subdir in subdirs:
        result = check_overfitting(subdir)
        if result == 1:
            overfitts += 1

    return overfitts

if __name__ == "__main__":
    logdir = r"I:\sf2026\Astrocytes_V2_run2"
    filecount = 0
    for subdir in Path(logdir).rglob("*"):
        if not subdir.is_dir():
            continue
        result = check_overfitting(subdir)
        print(f"{subdir}: {'Possible overfitting flagged.' if result == 1 else 'No overfitting flagged.'}")
        filecount += 1
        if filecount > 10:
            exit()