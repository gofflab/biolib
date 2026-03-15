"""Utilities for submitting and monitoring jobs on an IBM Platform LSF cluster.

Provides the LSFJob class for constructing, submitting, polling, killing, and
waiting on LSF batch jobs via the bsub/bjobs/bkill command-line tools.  Also
supports a 'local' pseudo-queue for running commands directly on the current
host without LSF.

Designed for use with Harvard's Odyssey LSF cluster but applicable to any
Platform LSF installation.
"""
import os
import re
import subprocess
import sys
import time

# from misc import pp  # rasmus library removed - not Python 3.12 compatible

#Constants
lsf_mem = 32
lsf_default_queue = "normal_parallel" # normal_parallel  since it has less users

#######################
#Error Handling
#######################
class LSFError(Exception):
	"""Exception raised for LSF-related errors.

	Attributes:
		value: String or object describing the error condition.
	"""
	def __init__(self,value):
		"""Initialises an LSFError with an error value.

		Args:
			value: A string or object describing the LSF error.
		"""
		self.value = value

	def __str__(self):
		"""Returns a string representation of the error value."""
		return repr(self.value)

#################
#Base Class
#################
class LSFJob(object):
	"""Represents a single LSF batch job with lifecycle management.

	Constructs the bsub command string, submits the job to LSF (or runs it
	locally), and provides methods to poll job status, wait for completion,
	and kill the job.

	Attributes:
		cmd_str: The shell command to execute.
		queue: LSF queue name (or 'local' for local execution).
		outfile: Path to the stdout capture file.
		errfile: Path to the stderr capture file.
		job_name: Optional LSF job name.
		group: Optional LSF job group.
		job_mem: Memory requirement in GB (capped at lsf_mem global).
		submit_flag: True after the job has been submitted.
		complete: True after the job has finished.
		status: Current job status string (e.g. 'PEND', 'RUN', 'DONE').
		jobID: LSF job ID integer (-999 before submission).
		submit_time: Submission timestamp from bjobs.
		exec_host: Host on which the job is/was running.
		submit_host: Host from which the job was submitted.
		bsub_str: List of tokens forming the complete bsub command.
	"""

	def __init__(self,cmd_str,job_name=None,job_group=None,blocking=False,outfilename=None,errfilename=None,queue_name=None,job_mem=None,job_cores=1,notify=None):
		"""Creates an LSFJob instance and constructs the bsub command.

		Args:
			cmd_str: The shell command string to submit as an LSF job.
			job_name: Optional LSF job name passed to bsub -J.
			job_group: Optional LSF job group passed to bsub -g.
			blocking: If True, add -K flag to bsub to block until job
				completes.  Avoid on Odyssey LSF (limiting resource).
			outfilename: Path for stdout redirection.  If None, a temporary
				file in 'tmp/' is created.
			errfilename: Path for stderr redirection.  If None, a temporary
				file in 'tmp/' is created.
			queue_name: LSF queue name.  Defaults to lsf_default_queue.
				Use 'local' to run without LSF.
			job_mem: Memory requirement in GB.  Capped at the module-level
				lsf_mem constant.
			job_cores: Number of cores requested (stored but not currently
				used in the bsub command).
			notify: If truthy, add -N flag to bsub to send email notification
				on job completion.
		"""
		self.cmd_str = cmd_str

		global lsf_default_queue
		if queue_name == None:
			self.queue = lsf_default_queue
		else:
			self.queue = queue_name

		if outfilename == None:
			self.outfile = tmp_name("bsub_out_")
		else:
			self.outfile = outfilename
		if errfilename == None:
			self.errfile = tmp_name("bsub_err_")
		else:
			self.errfile = errfilename

		self.job_name = job_name
		self.group = job_group
		self.job_mem = job_mem
		self.submit_flag = False
		self.complete = False
		self.status = 'NOT SUBMITTED'
		self.jobID= -999

		self.submit_time = ""
		self.exec_host = ""
		self.submit_host = ""

		bsub_str = ["bsub"]

		if notify:
			bsub_str.extend(["-N"])

		bsub_str.extend(["-q", self.queue])

		if self.job_name != None:
			bsub_str.extend(["-J", self.job_name])

		if self.group != None:
			bsub_str.extend(['-g', self.group])

		if blocking != False:
			bsub_str.extend(["-K"])

		global lsf_mem
		if job_mem != None and lsf_mem != None:
			self.job_mem = min(self.job_mem, lsf_mem)
			bsub_str.extend(["-R rusage[mem=%d]" % self.job_mem])

		bsub_str.extend(["-R span[hosts=1]"])

		bsub_str.extend(["-oo", self.outfile])
		bsub_str.extend(["-eo", self.errfile])
		bsub_str.extend(["%s" % self.cmd_str])

		self.bsub_str = bsub_str

		#Handle if queue == "local"
		if self.queue == "local":
			local_str = [""]
			local_str.extend([">", self.outfile])
			local_str.extend(["2>", self.errfile])

			#TODO: Add self.cmd_str to bsub_str so command actually gets run.
			self.bsub_str = local_str
			self.bsub_str.insert(0,self.cmd_str)

	def __repr__(self):
		"""Returns a verbose string representation including all attributes."""
		return "Instance of class LSF Job:\n\t%s\n\tSubmitted: %s\n\t Complete: %s\n" % (self.cmd_str,self.submit_flag,self.complete) + str(self.__dict__)

	def __str__(self):
		"""Returns the complete bsub command as a space-joined string."""
		return " ".join(self.bsub_str)

	def submit(self): # wait pend
		"""Submits the job to LSF (or runs it locally) and waits for it to enter a stable state.

		For LSF jobs, uses subprocess.Popen to call bsub, retrieves the job ID,
		and polls until the status transitions out of 'SUBMITTED'.  For local
		jobs, launches the process and returns immediately.

		Returns:
			0 on successful submission (or 0 for local job launch).

		Raises:
			LSFError: If the bsub command returns a non-zero exit code.
		"""
		if self.submit_flag == True:
			print("Job already submitted", file=sys.stderr)
			return 0# what do you return here?

		self.submit_proc = subprocess.Popen(self.bsub_str,shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

		#Handle local jobs
		if self.queue == "local":
			self.submit_flag = True
			self.status = 'RUN'
			self.submit
			self.jobID = self.submit_proc.pid
			print("Job running locally with pid %d" % self.jobID, file=sys.stderr)
			return 0

		#Handle queued jobs
		if self.submit_proc.wait() != 0:
			raise LSFError("Could not submit to LSF. Error %d" % self.submit_proc.poll())
		else:
			self.submit_flag = True
			self.status = 'SUBMITTED'
			self.submit_status = self.submit_proc.stdout.read().rstrip()
			self.getJobId()
			#Wait until job switched from submitted to pend/run
			while self.status in ['SUBMITTED'] :
				try:
					self.poll()
				except Exception as e:
					print('Exception poll error: %s\n' % e, file=sys.stderr)

		print(self.submit_status, file=sys.stderr)
		return self.submit_proc.wait()

	def poll(self):
		"""This will poll using bjobs for the specific jobID for a given instance of LSFJob"""
		if not self.submit_flag:
			return "Job not yet submitted"
		elif self.complete:
			return "Job completed"
		else:
			#Handle local jobs
			if self.queue == "local":
				if self.submit_proc.poll() == 0:
					self.complete = True
					self.status = 'DONE'
					return self.status
				if self.submit_proc.poll() == None:
					return self.status
				if self.submit_proc.poll() > 0 or self.submit_proc.poll() < 0:
					raise LSFError("Problem with local job %d. Error %d" % (self.jobID,int(self.submit_proc.poll())))
				return
			tmp = subprocess.Popen('bjobs -a -w %d' % self.jobID,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			tmp_err = tmp.stderr.read().rstrip()
			notfoundpat = re.compile(r"Job \<[0-9]+\> is not found")
			failedpat = "Failed in an LSF library call"

			#wait until the bjobs query returns  (not until the job itself is finished)
			while tmp.wait() > 0:
				if tmp_err.count(failedpat) > 0:
					print(tmp_err, file=sys.stderr)
					time.sleep(20)
					tmp = subprocess.Popen('bjobs -w %d' % self.jobID,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
					tmp_err = tmp.stderr.read().rstrip()
					continue
				else:
					raise LSFError(tmp_err) # never caught
			if tmp.wait() == 0:
				#If job is not found and was previously pending/running assume it is finished
				if notfoundpat.match(tmp_err):
					if self.status in ['RUN','PEND']:
						self.status = 'DONE'
						self.complete = True
						return self.status
					else: # was never run
						print("waited, job did not run " + tmp_err, file=sys.stderr)
						return tmp_err
				#else: job still exists, update its status
				tmp_lines = [x.rstrip() for x in tmp.stdout.readlines()]
				keys,values = [x.split() for x in tmp_lines]
				tmpDict = dict(zip(keys,values))
				#pp(tmpDict)
				self.status = tmpDict['STAT']
				self.submit_time = tmpDict['SUBMIT_TIME']
				self.exec_host = tmpDict['EXEC_HOST']
				self.submit_host = tmpDict['FROM_HOST']
				return self.status
			else:
				#Should not reach this line... CONSIDER erasing and doing while tmp.wait!=0
				raise LSFError("Problem with bjobs polling. Error %s" % tmp_err)

	def getJobId(self):
		"""Parses the LSF job ID from the bsub submission output.

		Extracts the integer job ID from the '<JOBID>' pattern in
		self.submit_status and stores it in self.jobID.  Prints a message to
		stdout if the job has not been submitted yet.
		"""
		if self.submit_flag:
			jobID_search = re.search(r"\<[0-9]+\>",self.submit_status)
			self.jobID = int(jobID_search.group().strip("><"))
			return
		else:
			print("Job not yet submitted.")
			return

	def kill(self):
		"""Kills the LSF job using bkill.

		Does nothing if the job has not been submitted or has no valid job ID.
		Loops until bkill returns 0, retrying if necessary.  On success, resets
		status to 'NOT SUBMITTED' and clears submit_flag and complete.
		"""
		#Added this to fix cases were kill fails because there is no job id
		if self.status in ['NOT SUBMITTED'] or self.jobID== -999 :
			self.status = 'NOT SUBMITTED'
			return
		tmp = subprocess.Popen('bkill %d' % self.jobID,shell=True) #This is not working : we kill but the tmp wait doesn't return 0
		while tmp.wait() > 0:
			time.sleep(3)
			if tmp.wait()< 0: #if were not able to kill , try again
				tmp = subprocess.Popen('bkill %d' % self.jobID,shell=True) #This is not working : we kill but the tmp wait doesn't return 0
		if tmp.wait() == 0:
			self.status = 'KILLED'
			self.submit_flag = False
			self.complete = False
			self.status = 'NOT SUBMITTED'
		return

	def wait(self):
		"""Blocks until the LSF job reaches a terminal state.

		Polls the job status every 30 seconds until status is no longer
		'SUBMITTED', 'PEND', 'RUN', or 'SUSP'.  Prints a warning to stderr
		if the job is suspended.  Sets status to 'DONE' and complete to True
		on exit.
		"""
		self.poll()
		if not self.submit_flag:
			print("Job not yet submitted")
			return
		while self.status in['SUBMITTED','PEND','RUN','SUSP']:
			time.sleep(30)
			self.poll()
			if self.status in ['SUSP']:
				print('SUSPENDED : %d \n' % self.jobID, file=sys.stderr)
		self.status = 'DONE'
		self.complete = True
		return


##############
#Helper functions
##############
def tmp_name(prefix):
	"""Generates a unique temporary file path inside a local 'tmp/' directory.

	Creates the 'tmp/' directory in the current working directory if it does
	not already exist, then returns a path of the form
	'tmp/<prefix><random_suffix>'.

	Args:
		prefix: String prefix for the temporary file name.

	Returns:
		A string file path for a temporary file that does not yet exist.
	"""
	import tempfile
	tmp_root = "tmp/"
	if os.path.exists(tmp_root):
		pass
	else:
		os.mkdir(tmp_root)
	return tmp_root + prefix + os.path.basename(tempfile.mktemp())
