"""Preserve soot-dem's original/custom targets when adding the fluid component."""
import json
from pathlib import Path
import subprocess
import sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from install_into_soot_dem import install,INCLUDE
ROOT=Path(__file__).resolve().parents[1]
ORIGINAL='''cmake_minimum_required(VERSION 3.20)
project(soot_afm)
set(MY_CUSTOM_SETTING preserved)
add_executable(restructuring src/restructuring.cpp)
add_executable(aggregation src/aggregation.cpp)
add_executable(my_custom_target src/custom.cpp)
'''
class MergeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.base=Path(self.tmp.name)
        self.repo=self.base/'soot';self.repo.mkdir();(self.repo/'src').mkdir()
        for f in ('aggregate.h','restructuring.cpp','aggregation.cpp','custom.cpp'):(self.repo/'src'/f).write_text('// original source\n')
        self.payload=self.base/'payload';(self.payload/'cmake').mkdir(parents=True)
        (self.payload/'cmake/SootCapillaryCoupling.cmake').write_text('# module\n')
        (self.payload/'cmake/CapillaryStandalone.cmake').write_text('# component build\n')
        (self.payload/'cmake/capillary_payload.json').write_text(json.dumps([{'source':'cmake/CapillaryStandalone.cmake','destination':'CMakeLists.txt'}]))
    def tearDown(self):self.tmp.cleanup()
    def git(self,*args):return subprocess.run(['git','-C',str(self.repo),*args],capture_output=True,text=True,check=True)
    def test_preserves_working_custom_targets_and_is_idempotent(self):
        file=self.repo/'CMakeLists.txt';file.write_text(ORIGINAL)
        result=install(self.repo,self.payload)
        self.assertEqual(result['original_cmake_source'],'working tree')
        self.assertTrue(file.read_text().startswith(ORIGINAL))
        self.assertEqual(file.read_text().count(INCLUDE),1)
        self.assertEqual((self.repo/'CMakeLists.txt.before-capillary-merge').read_text(),ORIGINAL)
        self.assertEqual((self.repo/'capillary/CMakeLists.txt').read_text(),'# component build\n')
        second=install(self.repo,self.payload);self.assertIsNone(second['backup'])
        self.assertEqual(file.read_text().count(INCLUDE),1)
        self.assertEqual((self.repo/'src/custom.cpp').read_text(),'// original source\n')
    def test_recovers_original_git_head_after_accidental_replacement(self):
        file=self.repo/'CMakeLists.txt';file.write_text(ORIGINAL)
        self.git('init');self.git('add','.');self.git('-c','user.name=Regression','-c','user.email=regression@example.invalid','commit','-m','Original build')
        standalone='project(sphinxsys_capillary_validation)\nadd_executable(pf_dem standalone.cpp)\n';file.write_text(standalone)
        result=install(self.repo,self.payload)
        self.assertEqual(result['original_cmake_source'],'Git HEAD')
        self.assertTrue(file.read_text().startswith(ORIGINAL))
        self.assertEqual((self.repo/'CMakeLists.txt.before-capillary-merge').read_text(),standalone)
    def test_recovers_earlier_original_when_replacement_was_committed(self):
        file=self.repo/'CMakeLists.txt';file.write_text(ORIGINAL)
        self.git('init');self.git('add','.');self.git('-c','user.name=Regression','-c','user.email=regression@example.invalid','commit','-m','Original build')
        file.write_text('project(sphinxsys_capillary_validation)\n')
        self.git('add','.');self.git('-c','user.name=Regression','-c','user.email=regression@example.invalid','commit','-m','Accidental replacement')
        result=install(self.repo,self.payload)
        self.assertTrue(result['original_cmake_source'].startswith('Git '))
        self.assertTrue(file.read_text().startswith(ORIGINAL))
    def test_does_not_replace_unknown_project(self):
        file=self.repo/'CMakeLists.txt';file.write_text('project(unrelated)\n')
        with self.assertRaisesRegex(ValueError,'refusing'):install(self.repo,self.payload)
        self.assertEqual(file.read_text(),'project(unrelated)\n')
        self.assertFalse((self.repo/'capillary').exists())
if __name__=='__main__':unittest.main()
