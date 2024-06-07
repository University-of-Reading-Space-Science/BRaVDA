import pytest
import startBravda
from datetime import datetime
import numpy as np
import os


def test_solar_min_run():
    """
    Function to run BRaVDA for the solar minimum period (CR2080) to test if the new BRaVDA produces the same output
    as the old BRaVDA.
    :return:
    """

    # finding the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # setting up the variables to run the BRaVDA function
    # forecast date - end date of the CR
    fcDate = datetime(2009, 3, 9)
    # specific configuration file for this run
    config_file = os.path.join(current_dir, 'BRaVDA_test', 'solar_min_config.dat')
    precondState = True
    # output directory
    output_dir = os.path.join(current_dir, 'BRaVDA_test', 'Solar_min_test')

    # running BRaVDA
    startBravda.bravdafunction(fcDate, config_file, obsToAssim=["STERA", "STERB", "OMNI"], usecustomens=False,
                               runoutputdir=output_dir, plottimeseries=True, corona='MAS', precondState=precondState)

    return


def test_solar_min_prior():
    """
    Function to test if the new BRaVDA prior is the same as the old BRaVDA prior for CR2080.
    :return:
    """

    # finding the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # loading in the reference prior
    prior_ref_file = os.path.join(current_dir, 'BRaVDA_test', 'Solar_min_test_reference', 'prior',
                                  'prior_MJDstart54871.txt')
    prior_ref = np.loadtxt(prior_ref_file)

    # loading in the new prior
    prior_new_file = os.path.join(current_dir, 'BRaVDA_test', 'Solar_min_test', 'prior', 'prior_MJDstart54871.txt')
    prior_new = np.loadtxt(prior_new_file)

    # testing to see if the priors are the same within some tolerance
    assert np.allclose(prior_ref, prior_new)

    return


def test_solar_min_post():
    """
    Function to test if the new BRaVDA posterior is the same as the old BRaVDA posterior for CR2080
    :return:
    """

    # finding the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # loading in the reference posterior
    post_ref_file = os.path.join(current_dir, 'BRaVDA_test', 'Solar_min_test_reference', 'posterior',
                                 'posterior_MJDstart54871.txt')
    post_ref = np.loadtxt(post_ref_file)

    # loading in the new posterior
    post_new_file = os.path.join(current_dir, 'BRaVDA_test', 'Solar_min_test', 'posterior',
                                 'posterior_MJDstart54871.txt')
    post_new = np.loadtxt(post_new_file)

    # testing to see if the priors are the same within some tolerance
    assert np.allclose(post_ref, post_new)

    return


def test_solar_max_run():
    """
    Function to run BRaVDA for the solar minimum period (CR2139) to test if the new BRaVDA produces the same output
    as the old BRaVDA.
    :return:
    """

    # finding the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # setting up the variables to run the BRaVDA function
    # forecast date - end date of the CR
    fcDate = datetime(2013, 8, 4)
    # specific configuration file for this run
    config_file = os.path.join(current_dir, 'BRaVDA_test', 'solar_max_config.dat')
    precondState = True
    # output directory
    output_dir = os.path.join(current_dir, 'BRaVDA_test', 'Solar_max_test')

    # running BRaVDA
    startBravda.bravdafunction(fcDate, config_file, obsToAssim=["STERA", "STERB", "OMNI"], usecustomens=False,
                               runoutputdir=output_dir, plottimeseries=True, corona='MAS', precondState=precondState)

    return


def test_solar_max_prior():
    """
    Function to test if the new BRaVDA prior is the same as the old BRaVDA prior for CR2080.
    :return:
    """

    # finding the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # loading in the reference prior
    prior_ref_file = os.path.join(current_dir, 'BRaVDA_test', 'Solar_max_test_reference', 'prior',
                                  'prior_MJDstart56480.txt')
    prior_ref = np.loadtxt(prior_ref_file)

    # loading in the new prior
    prior_new_file = os.path.join(current_dir, 'BRaVDA_test', 'Solar_max_test', 'prior', 'prior_MJDstart56480.txt')
    prior_new = np.loadtxt(prior_new_file)

    # testing to see if the priors are the same within some tolerance
    assert np.allclose(prior_ref, prior_new)

    return


def test_solar_max_post():
    """
    Function to test if the new BRaVDA posterior is the same as the old BRaVDA posterior for CR2080
    :return:
    """

    # finding the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # loading in the reference posterior
    post_ref_file = os.path.join(current_dir, 'BRaVDA_test', 'Solar_max_test_reference', 'posterior',
                                 'posterior_MJDstart56480.txt')
    post_ref = np.loadtxt(post_ref_file)

    # loading in the new posterior
    post_new_file = os.path.join(current_dir, 'BRaVDA_test', 'Solar_max_test', 'posterior',
                                 'posterior_MJDstart56480.txt')
    post_new = np.loadtxt(post_new_file)

    # testing to see if the priors are the same within some tolerance
    assert np.allclose(post_ref, post_new)

    return

