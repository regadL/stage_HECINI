from re import search


def retrieveAllAccuracy (path_filin):
    """
    Retrieve accuracy one global LDA
    args: -> path filin
    return: -> accuracy train
            -> accuracy test
            -> number of descriptors
    """
    
    filin = open (path_filin, "r")
    list_lines = filin.readlines ()
    filin.close ()
    nb_descriptor = list_lines[2].split (" ")[1].strip().replace ("\"", "")
    list_accuracy = []
    for line_file in list_lines : 
        if search("accuracy", line_file) : 
            list_accuracy.append (line_file.split (":")[1].strip().replace ("\"", ""))
    
    return nb_descriptor, list_accuracy[0], list_accuracy[1], list_accuracy[2]


    
def retrieveQualityLOO (path_filin):
    """
    Retrieve quality prediction in leave one out
    args: -> path filin
    return: -> accuracy 
            -> Recall
            -> precision
            -> sensibility
            -> specificity
    """
    filin = open (path_filin, "r")
    list_lines = filin.readlines ()
    filin.close ()
    
    nb_descriptor = list_lines[2].split (" ")[1].strip().replace ("\"", "")
    for line_file in list_lines : 
        if search("accuracy", line_file) : 
            accuracy = line_file.split (":")[1].strip().replace ("\"", "")
        elif search ("precision", line_file) : 
            precision = line_file.split (":")[1].strip().replace ("\"", "")
        elif search ("recall", line_file) : 
            recall = line_file.split (":")[1].strip().replace ("\"", "")
        elif search ("sensibility", line_file): 
            sensibility = line_file.split (":")[1].strip().replace ("\"", "")
        elif search ("specificity", line_file):
            specificity = line_file.split (":")[1].strip().replace ("\"", "")
    
    return nb_descriptor, accuracy, precision, recall, sensibility, specificity
 
 
 
 
    
