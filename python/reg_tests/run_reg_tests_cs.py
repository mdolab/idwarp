from __future__ import print_function
# =============================================================================
# Standard Python modules                                           
# =============================================================================
import os, sys, argparse

# =============================================================================
# Extension modules
# =============================================================================
import mdo_regression_helper as reg

# define scripts to run:
module_name = 'idwarp_cs'

# Get the optional commandline arguments:
parser = argparse.ArgumentParser()
parser.add_argument("-mode",default='compare',choices=['train','compare'],
                    help='Train generates reference file. Compare runs test')
parser.add_argument("-diff_cmd",default='xxdiff',
                    help='Command to run for displaying diff. Default: xxdiff')
parser.add_argument("-nodiff", action='store_true', help='Suppress\
 displaying the comparison if not successful')
parser.add_argument("-mpiexec",default='mpirun',
                    help='Command to use for running mpi')
args = parser.parse_args()

mode = args.mode
diff_cmd = args.diff_cmd
nodiff = args.nodiff
mpiexec= args.mpiexec

if mode == 'train':
    try:
        os.remove('%s_reg.ref'%(module_name))
        os.remove('%s_par_reg.ref'%(module_name))
    except OSError:
        pass

    # Run script: Note that the reference values are actually the
    # *REAL* values not the CS ones. 
    os.system('%s -np 1 python solve_script_cs.py >> %s_reg.ref 2>&1'%(
        mpiexec, module_name))

    # Run script in parallel: Note that the reference values are actually the
    # *REAL* values not the CS ones. 
    os.system('%s -np 4 python solve_script_cs.py >> %s_par_reg.ref 2>&1'%(
        mpiexec, module_name))
            
    # If we're training, we done (no comparison)
    sys.exit(0)
else:
    try:
        os.remove('%s_reg'%(module_name))
        os.remove('%s_par_reg'%(module_name))
    except OSError:
        pass

    # Run script
    os.system('%s -np 1 python solve_script_cs.py complex >> %s_reg 2>&1'%(
        mpiexec, module_name))

    # Run script
    os.system('%s -np 4 python solve_script_cs.py complex >> %s_par_reg 2>&1'%(
        mpiexec, module_name))

    # Do the comparison (reference file must be first)
    res = reg.reg_file_comp('%s_reg.ref'%(module_name),'%s_reg'%(module_name))
    res_par = reg.reg_file_comp('%s_par_reg.ref'%(module_name),'%s_par_reg'%(module_name))

# Set the proper return codes for the script running this:
if res == 0: #reg.REG_FILES_MATCH
    print('%s: Success!'%(module_name))
elif res == 1: #reg.REG_FILES_DO_NOT_MATCH
    print('%s: Failure!'%(module_name))
    if not nodiff:
        os.system('%s %s_reg.ref %s_reg'%(diff_cmd, module_name, module_name))
elif res == -1: #reg.REG_ERROR
    print('%s: Error in regression. Missing files.'%(module_name))

# Set the proper return codes for the script running this:
if res_par == 0: #reg.REG_FILES_MATCH
    print('%s: Parallel Success!'%(module_name))
elif res_par == 1: #reg.REG_FILES_DO_NOT_MATCH
    print('%s: Parallel Failure!'%(module_name))
    if not nodiff:
        os.system('%s %s_par_reg.ref %s_par_reg'%(diff_cmd, module_name, module_name))
elif res_par == -1: #reg.REG_ERROR
    print('%s: Error in parallel regression. Missing files.'%(module_name))


# Exit with code from reg_file_comp:
res = max(res,res_par)
sys.exit(res)

