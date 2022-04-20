from pexpect import popen_spawn

def add_ud_proteases():
    add_enzyme = popen_spawn.PopenSpawn('rpg -a')
    add_enzyme.expect('Name of the new enzyme')
    add_enzyme.sendline('Asp-N-UD')
    #add_enzyme.expect('Create a cleaving rule')
    add_enzyme.sendline('c')
    #add_enzyme.expect('Write your cleaving rule')
    add_enzyme.sendline('(,D)')
    add_enzyme.sendline('q')
    #add_enzyme.expect('Add another enzyme')
    add_enzyme.sendline('n')