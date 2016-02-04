import os
import datetime

PROJECT_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')


def s_timestamp():
    return datetime.datetime.now().strftime("%H_%M_%S %d_%m_%y")