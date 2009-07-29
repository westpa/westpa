from we.environment import *
        
def update_current_segment(*args, **kwargs):
    updates = {}
    for arg in args:
        updates.update(arg)
    updates.update(kwargs)

    segment = get_current_segment(load_parent = True)
    get_data_manager().update_segment(segment, updates)


        
    
    
    