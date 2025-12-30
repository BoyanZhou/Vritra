"""
set up logger for other script
"""

import logging
import os


def setup_log_file(log_file_path, log_level=logging.INFO, log_format=None):
    """
    Set up a log file with the specified output path.
    Args:
        log_file_path (str): The full path to the log file.
        log_level (int, optional): Logging level (e.g., logging.DEBUG, logging.INFO). Defaults to logging.INFO.
        log_format (str, optional): Format of the log messages. Defaults to a standard format.
    Returns:
        logging.Logger: Configured logger object.
    """
    # Ensure the directory for the log file exists
    if len(os.path.dirname(log_file_path)) > 0:
        os.makedirs(os.path.dirname(log_file_path), exist_ok=True)
    # Set a default log format if none is provided
    if log_format is None:
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    # Configure the logger
    logging.basicConfig(
        filename=log_file_path,
        level=log_level,
        format=log_format,
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    # Return the root logger
    return logging.getLogger()
