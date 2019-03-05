"""LOG file generation with error and monitoring of actions"""

# global module
import time
from re import sub

# personal module
import pathDirectory


def initAction( actionName ):
    """Initialisation of log file
    args: action
    return: file open and date Start"""
    
    dateStart = time.strftime( '%X %x' )
    dateStartFileName = formatdate( dateStart )
    print dateStart, ": ", actionName
    begin = time.time()
    fileLog = open( pathDirectory.logDir() + str( dateStartFileName ) + str( actionName.replace ( " ", "_" ) ), "w" )
    fileLog.write( actionName + "\n" )
    fileLog.write( dateStart + "\n" )

    return begin, fileLog


def endAction( action, start_time, log_file ):
    """
    End log file
    arg: - name action
         - time start
         - log file
    out: close file
    """
    
    dateEnd = time.strftime( '%X %x' )
    end = time.time()
    time_execution = end - start_time
    hour = int( time_execution / 3600 )
    timeRest = time_execution - ( hour * 3600 )
    minute = int( timeRest / 60 )
    second = int( timeRest - ( minute * 60 ) )
    print "Time execution :", hour, "h", minute, "min", second, "s"
    log_file.write( "Time execution :" + str( hour ) + "h" + str( minute ) + "min" + str( second ) + "s" + "\n" )
    print dateEnd, ": End ", action, "\n"
    log_file.write( dateEnd + ": End " + action )
    log_file.close()


def formatdate( date ):
    """
    Format date in time module
    arg: date
    return: date formated
    """
    
    date = sub( "/", "-", date )
    date = sub( " ", "_", date )

    return date + "_"

