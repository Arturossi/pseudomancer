import pytest
import os
import shutil
from unittest.mock import patch, MagicMock
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from pseudomancer.dependencies import check_dependencies
from pseudomancer.clustering import cluster_proteins_mmseqs2
from pseudomancer.search import run_mmseqs2_search


class TestDependencies:
    """Test dependency checking functionality."""
    
    def test_check_dependencies_available(self):
        """Test when all dependencies are available."""
        with patch('shutil.which') as mock_which:
            mock_which.return_value = '/usr/bin/mmseqs'
            # Should not raise an exception
            check_dependencies()
    
    def test_check_dependencies_missing(self):
        """Test when dependencies are missing."""
        with patch('shutil.which') as mock_which:
            mock_which.return_value = None
            with pytest.raises(SystemExit):
                check_dependencies()


class TestProteinClustering:
    """Test protein clustering functionality."""
    
    @patch('subprocess.run')
    def test_cluster_proteins_mmseqs2_success(self, mock_run, sample_protein_file, temp_dir):
        """Test successful protein clustering."""
        mock_run.return_value = MagicMock()
        
        # Create expected output file
        expected_output = os.path.join(temp_dir, "test_proteins_clustered_rep_seq.fasta")
        with open(expected_output, 'w') as f:
            f.write(">cluster_1\nMKRLLAISLLLAVVTSLLAAPYVKA\n")
        
        result = cluster_proteins_mmseqs2(sample_protein_file, temp_dir)
        
        # Check that mmseqs was called with correct parameters
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert call_args[0] == "mmseqs"
        assert call_args[1] == "easy-cluster"
        assert "--min-seq-id" in call_args
        assert "0.99" in call_args
        
        assert result == expected_output
    
    @patch('subprocess.run')
    def test_cluster_proteins_mmseqs2_failure(self, mock_run, sample_protein_file, temp_dir):
        """Test failed protein clustering."""
        from subprocess import CalledProcessError
        mock_run.side_effect = CalledProcessError(1, 'mmseqs', stderr="Error message")
        
        with pytest.raises(SystemExit):
            cluster_proteins_mmseqs2(sample_protein_file, temp_dir)


class TestMmseqs2Search:
    """Test mmseqs2 search functionality."""
    
    @patch('subprocess.run')
    def test_run_mmseqs2_search_success(self, mock_run, sample_protein_file, sample_genome_file, temp_dir):
        """Test successful mmseqs2 search."""
        mock_run.return_value = MagicMock()
        
        # Create expected output file
        results_dir = os.path.join(temp_dir, "mmseqs_results")
        os.makedirs(results_dir, exist_ok=True)
        expected_output = os.path.join(results_dir, "search_results.txt")
        with open(expected_output, 'w') as f:
            f.write("query\tqstart\tqend\tqcov\ttarget\ttstart\ttend\ttframe\tevalue\tpident\talnlen\tqlen\n")
            f.write("test_protein_1\t1\t25\t100\ttest_contig_1\t1\t75\t1\t1e-10\t80.0\t25\t25\n")
        
        result = run_mmseqs2_search(sample_protein_file, sample_genome_file, temp_dir)
        
        # Check that mmseqs commands were called
        assert mock_run.call_count >= 4  # createdb (2x), search, convertalis
        
        # Check some of the calls
        calls = [call[0][0] for call in mock_run.call_args_list]
        assert any("createdb" in call for call in calls)
        assert any("search" in call for call in calls)
        assert any("convertalis" in call for call in calls)
        
        assert result == expected_output
    
    @patch('subprocess.run')
    def test_run_mmseqs2_search_failure(self, mock_run, sample_protein_file, sample_genome_file, temp_dir):
        """Test failed mmseqs2 search."""
        from subprocess import CalledProcessError
        mock_run.side_effect = CalledProcessError(1, 'mmseqs', stderr="Database creation failed")
        
        with pytest.raises(CalledProcessError):
            run_mmseqs2_search(sample_protein_file, sample_genome_file, temp_dir)


class TestFileHandling:
    """Test file handling functionality."""
    
    def test_output_directories_created(self, temp_dir):
        """Test that output directories are created properly."""
        from pseudomancer.search import run_mmseqs2_search
        
        with patch('subprocess.run'):
            # Create mock files
            protein_file = os.path.join(temp_dir, "proteins.faa")
            genome_file = os.path.join(temp_dir, "genome.fasta")
            with open(protein_file, 'w') as f:
                f.write(">test\nMKR\n")
            with open(genome_file, 'w') as f:
                f.write(">test\nATG\n")
            
            try:
                run_mmseqs2_search(protein_file, genome_file, temp_dir)
            except:
                pass  # We expect this to fail without real mmseqs2
            
            # Check directories were created
            assert os.path.exists(os.path.join(temp_dir, "mmseqs_dbs"))
            assert os.path.exists(os.path.join(temp_dir, "mmseqs_results"))