if True:
	if True: print("hello")

if 0<1 and True: print("hello")

my_list = [1, 2, 3, 4, 5]
list_id_before = id(my_list)
print("ID before clearing:", list_id_before)

# Clear the list using one of the methods
my_list.clear()  # or my_list[:] = [] or del my_list[:]

list_id_after = id(my_list)
print("ID after clearing:", list_id_after)

# Verify that the ID has not changed
print("Reference unchanged:", list_id_before == list_id_after)  # Output: True

