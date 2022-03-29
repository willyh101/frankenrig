from psychopy import visual, core

win = visual.Window(size=[1920,1080], screen=1)
message = visual.TextStim(win, text='hello')
message.autoDraw = True
win.flip()
core.wait(2.0)
message.text = 'world'
win.flip()
core.wait(2.0)