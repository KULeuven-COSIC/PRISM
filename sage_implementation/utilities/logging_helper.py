import logging
import sys

# https://stackoverflow.com/questions/2183233/how-to-add-a-custom-loglevel-to-pythons-logging-facility/35804945#35804945

def addLoggingLevel(levelName, levelNum, methodName=None):
    """
    Comprehensively adds a new logging level to the `logging` module and the
    currently configured logging class.

    `levelName` becomes an attribute of the `logging` module with the value
    `levelNum`. `methodName` becomes a convenience method for both `logging`
    itself and the class returned by `logging.getLoggerClass()` (usually just
    `logging.Logger`). If `methodName` is not specified, `levelName.lower()` is
    used.
    To avoid accidental clobberings of existing attributes, this method will
    raise an `AttributeError` if the level name is already an attribute of the
    `logging` module or if the method name is already present 

    Example
    -------
    >>> addLoggingLevel('TRACE', logging.DEBUG - 5)
    >>> logging.getLogger(__name__).setLevel("TRACE")
    >>> logging.getLogger(__name__).trace('that worked')
    >>> logging.trace('so did this')
    >>> logging.TRACE
    5

    """
    if not methodName:
        methodName = levelName.lower()

    if hasattr(logging, levelName):
       raise AttributeError('{} already defined in logging module'.format(levelName))
    if hasattr(logging, methodName):
       raise AttributeError('{} already defined in logging module'.format(methodName))
    if hasattr(logging.getLoggerClass(), methodName):
       raise AttributeError('{} already defined in logger class'.format(methodName))

    def logForLevel(self, message, *args, **kwargs):
        if self.isEnabledFor(levelNum):
            self._log(levelNum, message, args, **kwargs)
    def logToRoot(message, *args, **kwargs):
        logging.log(levelNum, message, *args, **kwargs)

    logging.addLevelName(levelNum, levelName)
    setattr(logging, levelName, levelNum)
    setattr(logging.getLoggerClass(), methodName, logForLevel)
    setattr(logging, methodName, logToRoot)

addLoggingLevel('DEBUG1', logging.DEBUG - 1)
addLoggingLevel('DEBUG2', logging.DEBUG - 2)
addLoggingLevel('DEBUG3', logging.DEBUG - 3)
addLoggingLevel('INFO1', logging.INFO + 1)
addLoggingLevel('INFO2', logging.INFO + 2)
addLoggingLevel('INFO3', logging.INFO + 3)
addLoggingLevel('TIMINGS', logging.WARNING + 1)
addLoggingLevel('PRINT', logging.WARNING + 2)
addLoggingLevel('BOLD', logging.WARNING + 3)

# https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output/56944256#56944256
class CustomFormatter(logging.Formatter):
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    bold="\x1b[1m"
    reset = "\x1b[0m"
    #format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"
    format = "%(levelname)s:%(name)s:%(message)s:(%(filename)s:%(lineno)d)"
    simple_format = '%(message)s'

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.DEBUG1: grey + format + reset,
        logging.DEBUG2: grey + format + reset,
        logging.DEBUG3: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.INFO1: grey + format + reset,
        logging.INFO2: grey + format + reset,
        logging.INFO3: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.TIMINGS: simple_format,
        logging.PRINT: simple_format,
        logging.BOLD: bold + simple_format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

ch = logging.StreamHandler(sys.stdout)
#ch.setLevel(logging.DEBUG3)
ch.setFormatter(CustomFormatter())
logging.root.addHandler(ch)

if __name__ == "__main__":
    logger = logging.getLogger("Logging test")
    logger.setLevel(logging.DEBUG3)
    logger.debug("Debug")
    logger.debug1("Debug1")
    logger.debug2("Debug2")
    logger.debug3("Debug3")
    logger.info("Info")
    logger.info1("Info1")
    logger.info2("Info2")
    logger.info3("Info3")
    logger.warning("Warnings")
    logger.timings("Timings")
    logger.print("Print")
    logger.bold("Bold")
    logger.error("Error")
    logger.critical("Critical")
