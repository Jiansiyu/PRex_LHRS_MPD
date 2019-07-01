'''
'''
from fpdf import FPDF

class pdfGenerator(object):
    def __init__(self,pdfname="result/detector.pdf"):
        
        self.pdfname=pdfname
        self.pdf=FPDF()

    def generatorpdf(self,detectorID=1):
        
        # form the filename that ready to load 
        for eventID in range(1,100):
            self.pdf.add_page()
            self.pdf.set_font("Arial", size=12)
            self.pdf.cell(200, 10, txt="GEM{0} Event{1}".format(detectorID,eventID), ln=1,align="C")
            self.pdf.set_font("Arial", size=12)
            self.pdf.image('result/Canvas_event{0}_TS0_x{1}.jpg'.format(eventID,detectorID),x=2,y=18,w=100)
            self.pdf.image('result/Canvas_event{0}_TS1_x{1}.jpg'.format(eventID,detectorID),x=2,y=108,w=100)
            self.pdf.image('result/Canvas_event{0}_TS2_x{1}.jpg'.format(eventID,detectorID),x=2,y=198,w=100)
            self.pdf.image('result/Canvas_event{0}_TS0_y{1}.jpg'.format(eventID,detectorID),x=100,y=18,w=100)
            self.pdf.image('result/Canvas_event{0}_TS1_y{1}.jpg'.format(eventID,detectorID),x=100,y=108,w=100)
            self.pdf.image('result/Canvas_event{0}_TS2_y{1}.jpg'.format(eventID,detectorID),x=100,y=198,w=100)

        if detectorID<=3:
            self.pdf.output('SBU_GEM{0}.pdf'.format(detectorID))
        else:
            self.pdf.output('UVa_GEM{0}.pdf'.format(detectorID))

    def generatorAll(self):
        # form the filename that ready to load 
        for eventID in range(1,100):
            for detectorID in range(1,7):
                print('working on event {}'.format(eventID))
                self.pdf.add_page()
                self.pdf.image('result/Canvas_event{0}_TS0_x{1}.jpg'.format(eventID,detectorID),x=2,y=18,w=100)
                self.pdf.image('result/Canvas_event{0}_TS1_x{1}.jpg'.format(eventID,detectorID),x=2,y=108,w=100)
                self.pdf.image('result/Canvas_event{0}_TS2_x{1}.jpg'.format(eventID,detectorID),x=2,y=198,w=100)
                self.pdf.image('result/Canvas_event{0}_TS0_y{1}.jpg'.format(eventID,detectorID),x=100,y=18,w=100)
                self.pdf.image('result/Canvas_event{0}_TS1_y{1}.jpg'.format(eventID,detectorID),x=100,y=108,w=100)
                self.pdf.image('result/Canvas_event{0}_TS2_y{1}.jpg'.format(eventID,detectorID),x=100,y=198,w=100)
                self.pdf.set_font("Arial", size=12)
                self.pdf.cell(200, 2, txt="GEM{0} Event{1}".format(detectorID,eventID), ln=10,align="C")
                self.pdf.set_text_color(255,0,0)
                self.pdf.cell(10, 5, txt="Time Sample 0", ln=1,align="C")
                self.pdf.cell(10, 170, txt="Time Sample 1", ln=1,align="C")
                self.pdf.cell(10, 10, txt="Time Sample 2", ln=1,align="C")
        print("write to file")
        self.pdf.output('SBU_GEM1_3_50event.pdf')

    def addImage(self,image_path):
        pass

    
    

if __name__ == "__main__":
    test=pdfGenerator()
    test.generatorAll()

