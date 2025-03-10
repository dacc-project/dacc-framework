#!/usr/bin/env sage
# -*- coding: utf-8 -*-
"""
DACC Utilities Module

This module provides standardized utilities for the DACC framework, ensuring consistent
behavior across all components. It handles curve validation, LMFDB loading, and
configuration management to prevent inconsistencies in the analysis pipeline.

This centralized approach helps maintain the integrity of the research by ensuring
that all components use the same properly formatted curve labels and loading mechanisms.
"""
import warnings

# Suppress the specific urllib3 OpenSSL warning
warnings.filterwarnings("ignore", category=Warning, module="urllib3")

import requests
import re
import json
import os
from sage.all import EllipticCurve

def is_valid_lmfdb_label(label):
    """
    Check if a curve label has the proper LMFDB format.
    
    LMFDB labels follow the pattern 'conductor.isogeny_class'
    Examples: '11.a1', '389.a1', '5077.a1'
    """
    # Pattern: digits, dot, letter(s) followed by digit(s)
    pattern = r'^\d+\.[a-zA-Z]+\d+$'
    return bool(re.match(pattern, label))

def format_suggested_label(label):
    """
    Try to suggest a proper LMFDB format if the label is malformed.
    
    This function looks for the first letter in the string and
    inserts a dot before it to separate conductor from isogeny class.
    """
    if '.' in label:
        return label  # Already has a dot
        
    # Find where the letter part starts
    match = re.search(r'[a-zA-Z]', label)
    if match:
        letter_pos = match.start()
        return label[:letter_pos] + '.' + label[letter_pos:]
    
    return None  # Can't suggest a format

def get_curve_from_lmfdb(label, api_base="http://127.0.0.1:37777/api", debug=False):
    """
    Load a curve from LMFDB with proper validation.
    
    Args:
        label: The LMFDB label of the curve
        api_base: Base URL for the LMFDB API
        debug: Whether to print debug information
        
    Returns:
        (curve, label, data) tuple on success, or (None, None, None) on failure
    """
    if not is_valid_lmfdb_label(label):
        suggestion = format_suggested_label(label)
        print(f"ERROR: Invalid LMFDB label format: '{label}'")
        if suggestion:
            print(f"Did you mean: '{suggestion}'?")
        print("LMFDB labels should be in the format: 'conductor.isogeny_class'")
        print("Example: '11.a1', '389.a1', '5077.a1'")
        return None, None, None
        
    # Proceed with API call
    url = f"{api_base}/ec_curvedata/?lmfdb_label={label}&_format=json"
    
    if debug:
        print(f"Using URL: {url}")
        
    try:
        response = requests.get(url)
        
        if debug:
            print(f"Status code: {response.status_code}")
            print(f"Content-Type: {response.headers.get('Content-Type')}")
            
        if response.status_code == 200:
            data = response.json()
            
            if not data or 'data' not in data or not data['data']:
                print(f"Curve '{label}' not found in LMFDB database")
                print("Please check that the curve exists and the LMFDB service is accessible")
                return None, None, None
                
            curve_data = data['data'][0]
            
            # Create curve from ainvs
            ainvs = curve_data.get('ainvs')
            if not isinstance(ainvs, list):
                print(f"Invalid ainvs format in LMFDB data: {ainvs}")
                return None, None, None
                
            E = EllipticCurve(ainvs)
            return E, label, curve_data
            
        else:
            print(f"LMFDB request failed with status code: {response.status_code}")
            return None, None, None
            
    except Exception as e:
        print(f"Error retrieving curve from LMFDB: {e}")
        return None, None, None

def load_dacc_config(config_path="dacc_config.json"):
    """
    Load the DACC configuration file.
    
    Args:
        config_path: Path to the configuration file
        
    Returns:
        Configuration dictionary or default values if loading fails
    """
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        print(f"Loaded configuration from {config_path}")
        return config
    except Exception as e:
        print(f"Warning: Could not load configuration from {config_path}: {e}")
        print("Using built-in default configuration")
        
        # Minimal default config with just a few crucial curves
        return {
            "test_curves": {
                "by_rank": {
                    "0": ["11.a1"],
                    "1": ["37.a1"],
                    "2": ["389.a1"],
                    "3": ["5077.a1"],
                    "4": ["234446.a1"]
                }
            }
        }
    
def get_test_curves_by_rank(rank, limit=None, config=None):
    """
    Get test curves for a specific rank from configuration.
    
    Args:
        rank: The rank as an integer or string
        limit: Maximum number of curves to return
        config: Optional configuration dictionary (loaded from file if None)
        
    Returns:
        List of curve labels for the specified rank
    """
    if config is None:
        config = load_dacc_config()
        
    # Convert rank to string for dictionary lookup
    rank_str = str(rank)
    
    try:
        curves = config["test_curves"]["by_rank"].get(rank_str, [])
        if limit is not None and limit > 0:
            return curves[:limit]
        return curves
    except (KeyError, TypeError):
        print(f"Warning: Could not find test curves for rank {rank}")
        return []
        
def get_test_family(family_name, config=None):
    """
    Get a test family from configuration.
    
    Args:
        family_name: Name of the test family
        config: Optional configuration dictionary (loaded from file if None)
        
    Returns:
        Dictionary containing the test family
    """
    if config is None:
        config = load_dacc_config()
        
    try:
        return config["test_curves"]["test_families"].get(family_name, {})
    except (KeyError, TypeError):
        print(f"Warning: Could not find test family '{family_name}'")
        return {}
        
def get_default_examples(config=None):
    """
    Get default examples from configuration.
    
    Args:
        config: Optional configuration dictionary (loaded from file if None)
    
    Returns:
        Dictionary containing default examples for documentation and testing
    """
    if config is None:
        config = load_dacc_config()
        
    try:
        return config["test_curves"]["default_examples"]
    except (KeyError, TypeError):
        print("Warning: Could not find default examples")
        return {
            "sha_examples": [
                {"curve": "11.a1", "l_value": 0.254, "bsd_prediction": 0.254, "analytic_sha": 1}
            ],
            "rank_examples": {
                "1": [{"curve": "37.a1", "regulator": 0.051}]
            }
        }

def validate_curve_labels_in_config(config=None):
    """
    Validate all curve labels in the configuration to ensure they have proper LMFDB format.
    
    Args:
        config: Optional configuration dictionary (loaded from file if None)
        
    Returns:
        Tuple of (is_valid, issues) where issues is a list of problematic labels
    """
    if config is None:
        config = load_dacc_config()
        
    issues = []
    
    # Helper function to check labels recursively
    def check_labels_in_dict(d, path=""):
        for key, value in d.items():
            current_path = f"{path}.{key}" if path else key
            
            if isinstance(value, dict):
                check_labels_in_dict(value, current_path)
            elif isinstance(value, list):
                for i, item in enumerate(value):
                    if isinstance(item, dict):
                        check_labels_in_dict(item, f"{current_path}[{i}]")
                    elif isinstance(item, str) and re.search(r'\d+[a-zA-Z]+\d+', item):
                        # This looks like it might be a curve label
                        if not is_valid_lmfdb_label(item):
                            suggestion = format_suggested_label(item)
                            issues.append({
                                "path": f"{current_path}[{i}]",
                                "label": item,
                                "suggested": suggestion
                            })
            elif isinstance(value, str) and re.search(r'\d+[a-zA-Z]+\d+', value):
                # This looks like it might be a curve label
                if not is_valid_lmfdb_label(value):
                    suggestion = format_suggested_label(value)
                    issues.append({
                        "path": current_path,
                        "label": value,
                        "suggested": suggestion
                    })
    
    try:
        # Start the recursive check from the root
        check_labels_in_dict(config)
    except Exception as e:
        print(f"Error validating curve labels: {e}")
        
    return len(issues) == 0, issues