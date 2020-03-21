package abc;

import java.io.File;
import java.io.IOException;

import jxl.Workbook;
import jxl.write.Label;
import jxl.write.Number;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;
import jxl.write.WriteException;

public class ExcelIO {
    public static void main(String[] args) {
        ExcelIO excelIO = new ExcelIO();
        try {
            excelIO.JxlWriteDemo();
        } catch (IOException | WriteException e) {
            e.printStackTrace();
        }
    }


    private void JxlWriteDemo() throws IOException, WriteException {
        File xlsFile = new File("jxl.xls");
        // 创建一个工作簿
        WritableWorkbook workbook = Workbook.createWorkbook(xlsFile);
        // 创建一个工作表
        WritableSheet sheet = workbook.createSheet("sheet1", 0);

        for (int row = 0; row < 5; row++) {
            for (int col = 0; col < 5; col++) {
                // 向工作表中写入字符串
                Label label = new Label(col, row, "ss");
                sheet.addCell(label);
            }
            for (int col = 5; col < 10; col++) {
                // 写入数字
                Number number = new Number(col, row, 5);
                sheet.addCell(number);
            }
        }
        workbook.write();
        workbook.close();
    }
}
