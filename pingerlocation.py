import matlab.engine
import time
eng = matlab.engine.start_matlab()
# Call to pinger location algorithm
# First input is submarine depth, second input is required pinger frequency in kHz
location = eng.test(35)
print(location)

