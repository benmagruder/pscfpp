"""! Module for parsing state files produced a sweep. """

from pscfpp.param import *
from pscfpp.thermo import *
from pscfpp.state import *
import os

##
#   Container for data in state files produced by a PSCF sweep.
#
#   A Sweep is a container for all the data contained in the PSCF state
#   files created by a sweep. The data contained in each state file is 
#   stored in an instance of class State (full name pscfpp.state.State). 
#   Individual state objects within a Sweep can be accessed using a 
#   square-bracket syntax, like elements of a python list. Methods
#   summary and summaryString return summary reports containing values 
#   of selected variables.
#
#   **Construction:**
# 
#      The constructor for class Sweep parses and stores the contents 
#      of all of the the state files (which end with '.dat') produced 
#      by a PSCF sweep. The constructor takes a single argument which 
#      is a string that can be either the name of the parameter file
#      used to run the sweep, or the name of the directory in which
#      the resulting state files are located.
#
#      Examples: 
#
#      1. By passing in the parameter file name 'param':
#         \code
#            from pscfpp.sweep import *
#            s = Sweep('param')
#         \endcode
#
#      2. By passing in a directory of the folder that contains all 
#         the 'state files', with name 'out/':
#         \code
#            from pscfpp.sweep import *
#            s = Sweep('out/')
#         \endcode
#
#   **Accessing elements:**
#      
#      A Sweep object can be treated as a list of State objects, 
#      each of which contains the contents of a corresponding state 
#      file. Each state object can be accessed by a corresponding 
#      index with square brackets, like an element of a list, such
#      s[0] or s[1]. 
#
#      All elements and properties in each such State object can be
#      accessed using the dot notation for attributes of a State object,
#      such as s[0].param.Mixture.nMonomer or s[1].thermo.fHelmholtz.
#      See the documentation for class pscfpp.state.State for details.
#
#
#  **Summary reports:** 
#
#     See documentation for methods summary and summaryString.
#
class Sweep:

   ##
   # Constructor.
   #
   # The input parameter d is a string that can be either the name of the
   # directory that contains all of the data files or the name of the
   # parameter file used to perform the sweep. If the latter, the name
   # of the directory is extracted from the parameter file. 
   #
   # \param d  string that is the name of a directory or parameter file.
   #
   def __init__(self, d):
      old_cwd = os.getcwd()

      if os.path.isdir(d) == True:
         os.chdir(d)
      elif os.path.isfile(d) == True:
         p = Composite(d)
         for x,y in p.children.items():
            if x.rfind('Sweep') != -1:
               nd = old_cwd + '/' + y.baseFileName
               break
         try:
            os.chdir(nd)
         except:
            raise Exception('No Sweep block in the Param file')
      else:
         raise Exception('Not valid argument')

      filenames = []
      for file in os.listdir():
         if file.endswith('.dat'):
            filenames.append(file)
      filenames.sort()

      self.sweep = []
      for i in range(0, len(filenames)):
         s = State(filenames[i])
         self.sweep.append(s)

      os.chdir(old_cwd)

   ##
   # Make a summary report containing values for selected variables.
   #
   # This function constructs a data structure containing the values 
   # of a user-specified set of variables at each state in the sweep. 
   #
   # Specifying a list of variable names:
   #
   # The input parameter l is a list of strings that specify the 
   # variables of interest. Each element in this list is a string that 
   # gives the the name of an attribute of a single State object, starting 
   # with either 'param' or 'thermo', using a dot notation. For example
   # "param.Mixture.nMonomer" specifies the input parameter nMonomer of 
   # the Mixture block of the parameter section of the state file, while
   # "thermo.fHelmoltz" is the free energy per monomer reference volume.
   #
   # Return value:
   # 
   # The quantity returned by this function is a list in which each 
   # element is itself a list of values of values of the specified 
   # variables at a single state point.  If the return value is assigned 
   # to a variable named "summary" then summary[2] is a list of values 
   # of the chosen variables at state 2. 
   #
   # The state index (e.g., 2) may optionally be added as the first 
   # element of this list of variable values, by setting the parameter 
   # index to True. By default index = False, so no index numbers are
   # included when this parameter is absent. 
   #
   # Example:  
   #
   # If s is a Sweep object representing a mixture that has been 
   # subjected to a sweep that changes composition, then the command
   # \code
   #   vals = s.summary(['param.Mixture.Polymer[0].phi','thermo.pressure'])
   # \endcode
   # constructs list named vals of values for the volume fraction for
   # the first polymer species and the non-dimensionalized pressure of 
   # the system, for each state in the sweep.  In this example, the 
   # resulting quantity vals is a list of lists something like
   # \code
   #            [[0.5, 32.4415250701], [0.55, 30.6782301376], 
   #             [0.6, 28.8344576024], [0.65, 26.8019750534],
   #             [0.7, 24.5205812724]]
   # \endcode
   # In this example, the vlue of vals[2] is a list of two values, given
   # by [0.6, 28.8344576024], in which the first element (0.6) is the 
   # volume fraction of polymer 0 and the second element (28.8344) is 
   # the value of the pressure, both evaluated in state 2 (the third 
   # state in the sweep). 
   #
   # Example:
   #
   # The command 
   # \code
   #    s.summary(['param.Interaction.chi[0][1]'], index = True)
   # \endcode
   # applied to a sweep that changes the indicated chi parameter will
   # return an array that may look something like the following:
   # \code
   #   [[0, 12.0], [1, 13.0], [2, 14.0], [3, 15.0], [4, 16.0]]
   # \endcode
   # Here the first value for each state is the integer state index, and 
   # the second is the value of the chi parameter chi[0][1] for 
   # interactions between monomers of types 0 and 1. 
   #
   # \param l list of strings that contain variable names in dot notation
   # \param index if True, add the state index as the first variable
   #
   def summary(self, l, index = False):
      n = len(l)

      summary = []

      for i in range(0, len(self.sweep)):
         if index == True:
            s = [i]
         else:
            s = []
         for j in range(0, n):
            a = self.sweep[i]
            string = 'a.' + l[j]
            try:
               val = eval(string)
            except RecursionError:
               raise Exception('Wrong command or values do not exist')
            s.append(val)
         summary.append(s)

      return summary

   ##
   # Return a summary report as a formatted string suitable for printing.
   #
   # This method produced a formatted string containing the same type of
   # summary report that is generated by the summary method, giving values
   # of selected variables for every state in a sweep. The parameter l is
   # a list of variable name strings formatted in dot notation, exactly as
   # for the corresponding parameter of the Sweep.summary method. The 
   # format includes labels for the variable names and includes a state
   # index for each state.
   #
   # Example: If applied to a sweep with 5 state points, the command
   # \code
   #    report = s.summaryString(['param.Interaction.chi[0][1]',\
   #                                 'thermo.fHelmholtz']))
   #    print(report)
   # \endcode
   # yields an output that looks somethink like this. 
   # 
   # \code       
   #        step      chi[0][1]     fHelmholtz
   #           0  1.2000000e+01  1.9256750e+00
   #           1  1.3000000e+01  2.1102042e+00
   #           2  1.4000000e+01  2.2716872e+00
   #           3  1.5000000e+01  2.4158122e+00
   #           4  1.6000000e+01  2.5464487e+00
   # \endcode
   # The string returned by summaryString can also be written to a file.
   #
   # \param l list of strings that contain variable names in dot notation
   #
   def summaryString(self, l):
      n = len(l)

      summary = []

      for i in range(0, len(self.sweep)):
         s = [i]
         for j in range(0, n):
            a = self.sweep[i]
            string = 'a.' + l[j]
            try:
               val = eval(string)
            except RecursionError:
               raise Exception('Wrong command or values do not exist')
            s.append(val)
         summary.append(s)

      nl = [4]
      nameList = ['step']
      for i in range(0, n):
         index = l[i].rfind('.')
         name = l[i][index+1:]
         nameList.append(name)
         nl.append(len(name))

      valType = []
      for i in range(0, len(summary)):
         valType.append([])
         for j in range(0, len(summary[0])):
            valType[i].append(type(summary[i][j]))

      for i in range(0, len(summary)):
         for j in range(0, len(summary[0])):
            length = len(str(summary[i][j]))
            if (valType[i][j] == str) and (length > nl[j]):
               nl[j] = length
            if (valType[i][j] == float) and (13 > nl[j]):
               nl[j] = 13
            if (valType[i][j] == int) and (length > nl[j]):
               nl[j] = length


      summaryString = ' '
      for i in range(0, len(nameList)):
         stringFormat = '{:>' + str(nl[i]) + 's}'
         if i != len(nameList)-1:
            summaryString += stringFormat.format(nameList[i]) + '  '
         else:
            summaryString += stringFormat.format(nameList[i]) + '\n '

      for i in range(0, len(summary)):
         for j in range(0, len(summary[0])):
            if valType[i][j] == int:
               stringFormat = '{:>' + str(nl[j]) + '}'
               summaryString += stringFormat.format(summary[i][j])
            if valType[i][j] == float:
               stringFormat = '{:.7e}'
               val = stringFormat.format(summary[i][j])
               summaryString += val
            if valType[i][j] == str:
               stringFormat = '{:>' + str(nl[j]) + 's}'
               summaryString += stringFormat.format(summary[i][j])
            if j == len(summary[0])-1:
               summaryString += '\n '
            else:
               summaryString += '  '

      return summaryString

   ##
   # Get a specific State object, specified by integer index.
   #
   # \param key  integer state index.
   #
   def __getitem__(self, key):
      return self.sweep[key]

   ##
   # Get the number of states in a sweep.
   #
   def __len__(self):
      return len(self.sweep)

