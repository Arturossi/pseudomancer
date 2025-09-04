import pytest
import os
import sys
from unittest.mock import patch, MagicMock
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from pseudomancer.__main__ import get_args, main


class TestArgumentParsing:
    """Test command line argument parsing."""
    
    def test_get_args_valid_input(self):
        """Test parsing valid arguments."""
        test_args = [
            '--genome', '/path/to/genome.fasta',
            '--taxon', 'Mycobacterium',
            '--out_dir', '/path/to/output',
            '--evalue', '1e-6'
        ]
        
        with patch('sys.argv', ['pseudomancer'] + test_args):
            args = get_args()
            
            assert args.genome_file == '/path/to/genome.fasta'
            assert args.taxon == 'Mycobacterium'
            assert args.output_dir == '/path/to/output'
            assert args.e_value == 1e-6
    
    def test_get_args_missing_required(self):
        """Test parsing with missing required arguments."""
        test_args = ['--genome', '/path/to/genome.fasta']
        
        with patch('sys.argv', ['pseudomancer'] + test_args):
            with pytest.raises(SystemExit):
                get_args()
    
    def test_get_args_default_evalue(self):
        """Test that default e-value is set correctly."""
        test_args = [
            '--genome', '/path/to/genome.fasta',
            '--taxon', 'Mycobacterium',
            '--out_dir', '/path/to/output'
        ]
        
        with patch('sys.argv', ['pseudomancer'] + test_args):
            args = get_args()
            assert args.e_value == 1e-5


class TestMainFunction:
    """Test the main function execution."""
    
    @patch('pseudomancer.__main__.run_pseudomancer_pipeline')
    @patch('pseudomancer.__main__.check_dependencies')
    @patch('os.path.exists')
    @patch('os.makedirs')
    def test_main_execution(self, mock_makedirs, mock_exists, mock_check_deps, mock_pipeline):
        """Test main function calls all necessary components."""
        mock_exists.return_value = False
        mock_pipeline.return_value = "/path/to/results.txt"
        
        test_args = [
            'pseudomancer',
            '--genome', '/path/to/genome.fasta',
            '--taxon', 'Mycobacterium',
            '--out_dir', '/path/to/output',
            '--evalue', '1e-6'
        ]
        
        with patch('sys.argv', test_args):
            main()
        
        # Check that dependencies were checked
        mock_check_deps.assert_called_once()
        
        # Check that output directory was created
        mock_makedirs.assert_called_once_with('/path/to/output', exist_ok=True)
        
        # Check that pipeline was called with correct arguments
        mock_pipeline.assert_called_once_with(
            taxon='Mycobacterium',
            genome_file='/path/to/genome.fasta',
            output_dir='/path/to/output',
            evalue=1e-6
        )
    
    @patch('pseudomancer.__main__.run_pseudomancer_pipeline')
    @patch('pseudomancer.__main__.check_dependencies')
    @patch('os.path.exists')
    @patch('os.makedirs')
    def test_main_existing_output_dir(self, mock_makedirs, mock_exists, mock_check_deps, mock_pipeline):
        """Test main function when output directory already exists."""
        mock_exists.return_value = True
        mock_pipeline.return_value = "/path/to/results.txt"
        
        test_args = [
            'pseudomancer',
            '--genome', '/path/to/genome.fasta',
            '--taxon', 'Mycobacterium',
            '--out_dir', '/path/to/output'
        ]
        
        with patch('sys.argv', test_args):
            main()
        
        # Check that makedirs was NOT called since directory exists
        mock_makedirs.assert_not_called()