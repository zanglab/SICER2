import sys
import os
import logging
import time

name = "sicer"


def set_log_path(log_path):
    RUN_LEVEL_NUM = 25
    SETUP_LEVEL_NUM = 35
    logging.addLevelName(RUN_LEVEL_NUM, "RUN")
    logging.addLevelName(SETUP_LEVEL_NUM, "SETUP")

    def log_run(self, message, *args, **kwargs):
        if self.isEnabledFor(RUN_LEVEL_NUM):
            self._log(RUN_LEVEL_NUM, message, args, **kwargs)

    def log_setup(self, message, *args, **kwargs):
        if self.isEnabledFor(SETUP_LEVEL_NUM):
            self._log(SETUP_LEVEL_NUM, message, args, **kwargs)

    logging.Logger.run = log_run
    logging.Logger.setup = log_setup

    logger = logging.getLogger("SICER 2")

    if not logger.hasHandlers() or len(logger.handlers) == 0:
        logger.propagate = False
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s", datefmt="%H:%M:%S")
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)

        # log with timestamp
        log_file = os.path.join(log_path, f'execution_log_{time.strftime("%Y%m%d-%H%M%S")}.txt')
        if os.path.exists(log_file):
            os.remove(log_file)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)

        logger.addHandler(console_handler)
        logger.addHandler(file_handler)

    simple_logger = logging.getLogger("s_logger")
    simple_formatter = logging.Formatter("%(message)s")
    simple_handler = logging.StreamHandler(sys.stdout)
    simple_handler.setFormatter(simple_formatter)
    simple_logger.addHandler(simple_handler)
    simple_file_handler = logging.FileHandler(log_file)
    simple_file_handler.setFormatter(simple_formatter)
    simple_logger.addHandler(simple_file_handler)
    simple_logger.setLevel(logging.INFO)
    simple_logger.propagate = False

    return logger, simple_logger
